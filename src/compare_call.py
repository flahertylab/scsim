# Functions for comparing VCF variant calls files

import os
import vcf
import pandas as pd
import pdb

def read_true(refdir, NSC):
    """Reads BED files containing true SNV locations and returns DataFrame
    with columns {sample, loc, snv}."""

    true_df = pd.DataFrame(columns=['sample', 'loc', 'snv'])
    for i in range(1, NSC + 1):
        filepath = os.path.join(refdir, 'snv_sc%d.bed' % i)
        ds = pd.read_csv(filepath, sep='\t').iloc[:,2]
        df = pd.DataFrame({'sample': 'after_wga_sc%d' % i, 'loc': ds, 'snv': True})
        true_df = true_df.append(df, ignore_index=True, sort = True)

    return true_df

def read_monovar(callfile='calls/monovar/all.vcf'):
    """Read VCF file containing called SNV locations from Monovar and returns
    DataFrame with columns {sample, loc, monovar}."""

    monovar_df = pd.DataFrame(columns=['sample', 'loc', 'monovar'])
    monovar_vcf_reader = vcf.Reader(open(callfile, 'r'))
    for record in monovar_vcf_reader:
        for sample in monovar_vcf_reader.samples:
            call=record.genotype(sample)
            monovar_df = monovar_df.append({'sample': call.sample[:-4], 'loc': record.POS-1, 'monovar': call.called}, ignore_index=True)

    return monovar_df

def read_bcftools(vcf_list):

    bcftools_df = pd.DataFrame(columns=['sample', 'loc', 'bcftools'])
    for callfile in vcf_list:
        bcftools_vcf_reader = vcf.Reader(open(callfile, 'r'))
        for record in bcftools_vcf_reader:
            bcftools_df = bcftools_df.append({'sample': os.path.basename(callfile[:-4]), 'loc': record.POS - 1, 'bcftools': True},
                                             ignore_index=True)

    return bcftools_df

def read_snver(vcf_list):
    """Read SNVer VCF file into pandas dataframe."""

    snver_df = pd.DataFrame(columns=['sample', 'loc', 'snver'])
    for callfile in vcf_list:
        snver_vcf_reader = vcf.Reader(open(callfile, 'r'))
        for record in snver_vcf_reader:
            snver_df = snver_df.append({'sample': os.path.basename(callfile[:-11]), 'loc': record.POS - 1, 'snver': True}, ignore_index=True)

    return snver_df

def combine_dfs(df1, df2):
    """Returns outer join of two DataFrames on {sample, loc}."""
    
    #Ensure that the loc column is of type int64
    df1 = df1.astype(dtype= {"loc":"int64"})
    df2 = df2.astype(dtype= {"loc":"int64"})
    
    # Join the true and called locations
    combined_df = pd.merge(df1, df2,
                            on=['sample', 'loc'],
                            how='outer',
                            sort=True)
    combined_df = combined_df.astype(dtype= {"loc":"int64",
                                             "sample":"object"})
    combined_df.fillna(False, inplace=True)

    return combined_df

