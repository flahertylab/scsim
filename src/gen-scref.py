from Bio import SeqIO
import numpy as np
import scipy.stats as ss
import pandas as pd
import time
import sys, os
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

MAX_ITER = 20
RES_DIR = ''

def filter_source(source=None, start=10**5, target_size=10**6,
                  target_format='fasta'):
    """ Function to read and filter fasta file  """

    # convert fasta to character string with BioPython
    for seq in SeqIO.parse(source, target_format):
        chr_full = str(seq.seq).upper()  # convert all base pairs to upper case

    source_size = len(chr_full)

    if target_size > source_size:
        raise ValueError('Target size exceeds source size! \
        Please input a smaller target size.')

    else:
        end = start + target_size

        chr_filt = chr_full[start:end] # filter source to fit target size

        logger.debug('Source length = %d. Sample length = %d.'%(len(chr_full),
                                                                len(chr_filt)))

        return(chr_filt)


def write_snv_mtx(snv_mtx, snv_loc, e):
    """ Function to store SNV flags matrix in csv and LaTeX formats """

    (n_snv, n_sc) = np.shape(snv_mtx)

    # cast flags matrix as pandas dataframe
    snv_df = pd.DataFrame(snv_mtx)

    # rename dataframe rows and columns
    snv_df.columns = ["SC%d" % i for i in range(1,n_sc+1)]
    snv_df.index = snv_loc+1

    # write data frame to markdown
    snv_df.to_csv(path_or_buf = os.path.join(RES_DIR, 'exp%d_snv_flags.md'%(e+1)), sep='|', index=True, header=True, mode='w', float_format='%.2g')
    
    # export data frame to latex
    snv_df.to_latex(os.path.join(RES_DIR, 'exp%d_snv_flags.tex'%(e+1)))


    
def write_fasta(seq, desc):
    """ Utility function to write sequence in fasta format """
    
    # set outfile name
    out=desc+'.fa'
    
    # write header and sequence to file
    fasta_f = open(os.path.join(RES_DIR,out), "w")
    fasta_f.write(">" + str(desc) + "\n" + ''.join(seq) + "\n")
    fasta_f.close()    
    
def write_sc(sc_obj, seq_type='wga_gtype', out="output"):
    """ Utility function to write a sequence to two fasta files,
    one for each allele"""
    
    # extract sequence for each allele 
    N = len(sc_obj)
    
    seq_a1 = []
    seq_a2 = []
    
    # loop throuch loci
    for i in range(N):         
        # extract gnotype as list
        gtype = list(sc_obj[i][seq_type])
        seq_a1.append(gtype[0])
        seq_a2.append(gtype[1])
            
    # write to two fasta files, one for each allele
    write_fasta(seq=seq_a1, desc=out+'_a1')
    write_fasta(seq=seq_a2, desc=out+'_a2')

def apply_ado(src_gt):
    """ Utility function to apply ADO on a given genotype """
    return(2*src_gt[np.random.choice(2)])

def apply_fp(src_gt):
    """Utility function to apply fp on a given genotype """
    
    tmp = list(src_gt)
    fp_gt = None
    # assert homozygous genotype
    if tmp[0] == tmp[1]:
        if tmp[0] == 'A':
            fp_gt = 'T'+tmp[1]
        elif tmp[0] == 'T':
            fp_gt = 'A'+tmp[1]
        elif tmp[0] == 'G':
            fp_gt = 'C'+tmp[1]
        elif tmp[0] == 'C':
            fp_gt = 'G'+tmp[1]
        else:
            logger.info('Unrecognized genotype.')
    else:
        fp_gt = src_gt
        
    return(fp_gt)


def gen_snv_transition_matrix(ts_tv_p=0.71,
                              ts_hm_ht_p=0.5,
                              tv_l_p=0.5,
                              tv_hm_ht_p=0.8):
    """ Function to compute SNV matrix of transition probabilities,
    default rates follow SInC reference (Fig1.a) """

    G_ref = {'AA':0, 'CC':1, 'GG':2, 'TT':3}
    G_alt = {0:'AA',1:'AC',2:'AG',3:'AT',4:'CC',5:'CG',6:'CT',7:'GG',8:'GT',9:'TT'}

    # initialize matrices of transition probabilities
    T_transition = np.zeros((len(G_ref),len(G_alt)))
    T_transversion_r = np.zeros((len(G_ref),len(G_alt)))
    T_transversion_l = np.zeros((len(G_ref),len(G_alt)))
    
    # fill in transition probabilties with specified rates
    T_transition[0,7] = T_transition[1,9] = T_transition[2,0] = T_transition[3,4] = ts_hm_ht_p
    T_transition[0,2] = T_transition[1,6] = T_transition[2,2] = T_transition[3,6] = 1 - ts_hm_ht_p

    T_transversion_l[0,4] = T_transversion_l[1,0] = T_transversion_l[2,9] = T_transversion_l[3,7] = tv_hm_ht_p
    T_transversion_l[0,1] = T_transversion_l[1,1] = T_transversion_l[2,8] = T_transversion_l[3,8] = 1 - tv_hm_ht_p

    T_transversion_r[0,9] = T_transversion_r[1,7] = T_transversion_r[2,4] = T_transversion_r[3,0] = tv_hm_ht_p
    T_transversion_r[0,3] = T_transversion_r[1,5] = T_transversion_r[2,5] = T_transversion_r[3,3] = 1 - tv_hm_ht_p
    
    # aggregate intermediate matrices
    T_transversion = tv_l_p * T_transversion_l + (1 - tv_l_p) * T_transversion_r 
    T_snv = ts_tv_p * T_transition + (1 - ts_tv_p) * T_transversion
    
    # assert final matrix rows sum to 1
    if np.sum(T_snv) != len(G_ref):
        raise ValueError('Transition probabilities matrix rows do not sum to 1!')
    
    return(T_snv)
        
        
def sim_alt_gtype(ref, T_snv):
    """ Function to simulate alternate SNV genotype using a matrix of
    transition probabilities """
    
    G_ref = {'AA':0, 'CC':1, 'GG':2, 'TT':3}
    G_alt = {0:'AA', 1:'AC', 2:'AG', 3:'AT', 4:'CC', 5:'CG', 6:'CT', 7:'GG', 8:'GT', 9:'TT'}
    
    # get row of probabilities corresponding to reference input
    ref_p = T_snv[G_ref[ref],:]
    # sample alternate genotype from reference probability vector
    alt_sample = ss.multinomial.rvs(1, p=ref_p)
    alt_idx = np.where(alt_sample == 1)[0][0] # get sample index
    alt_gtype = G_alt[alt_idx] # convert index to genotype
    
    return(alt_gtype)
    

def simulate_snv(source, snv_loc, alt_gtype, ADO_p = 0.2, FP_p = 3.2e-5):
    """ Function to simulate all SNVs in given locations with specified
    rates """

    n_snv = len(snv_loc)
    N = len(source)

    sc_gtypes = [dict(loc=i,
                      ref_gtype=2*source[i],
                      isSNV=False,
                      alt_gtype=2*source[i],
                      isADO=False,
                      ado_gtype=2*source[i],
                      isFP=False,
                      fp_gtype=2*source[i],
                      wga_gtype=2*source[i])
                 for i in range(N)]

    # Apply SNVs
    for i in snv_loc:
        sc_gtypes[i]['isSNV'] = True
        sc_gtypes[i]['alt_gtype'] = alt_gtype[i]
        sc_gtypes[i]['ado_gtype'] = sc_gtypes[i]['alt_gtype']
        sc_gtypes[i]['fp_gtype'] = sc_gtypes[i]['ado_gtype']
        sc_gtypes[i]['wga_gtype'] = sc_gtypes[i]['fp_gtype']

    # Apply ADOs
    ado_loc = np.random.choice(N, int(ADO_p*N))
    for i in ado_loc:
        sc_gtypes[i]['isADO'] = True
        sc_gtypes[i]['ado_gtype'] = apply_ado(sc_gtypes[i]['alt_gtype'])
        sc_gtypes[i]['fp_gtype'] = sc_gtypes[i]['ado_gtype']
        sc_gtypes[i]['wga_gtype'] = sc_gtypes[i]['fp_gtype']

    # Apply FPs
    fp_loc = np.random.choice(N, int(FP_p*N))
    for i in fp_loc:
        sc_gtypes[i]['isFP'] = True
        sc_gtypes[i]['fp_gtype'] = apply_fp(sc_gtypes[i]['ado_gtype'])
        sc_gtypes[i]['wga_gtype'] = sc_gtypes[i]['fp_gtype']

    return(sc_gtypes)


def sim_snv_flags_mtx(n_snv, n_sc):
    """ Function to simulate SNV flags matrix """

    # initialize SNV flags matrix
    snv_mtx = np.zeros((n_snv,n_sc), dtype=bool)

    # prepare indices for three simulation scenarios
    row_idx = np.linspace(start=0, stop=n_snv, num=4).astype(int)
    col_idx = int(n_sc/2)

    # scenario (1): all single cells share one third of the SNV locations
    snv_mtx[:row_idx[1],] = True

    # scenario (2): half of the single cells share one third
    # of the SNV locations
    snv_mtx[row_idx[1]:row_idx[2], col_idx:] = True

    # scenario (3): simulate shared and singleton SNVs uniformly
    p = {} # initialize SNV simulation proportions
    for i in range(row_idx[2],row_idx[3]):
        # simulate and store proportion of shared SNVs across all
        # single cells in row i
        p[i] = ss.uniform.rvs(0,1)
        # simulate SNV locations (encoded as '1's) for all single cells
        # in row i, based on simulated proportion p_i
        snv_mtx[i] = ss.binom.rvs(1, p[i], size=n_sc) == 1

    return(snv_mtx, p)


def simulate_sc(source, n_snv, n_sc, snv_loc=None, out='output', ts_tv_p=0.71, ts_hm_ht_p=0.5, AGCT_AGTC_p = 0.5, tv_hm_ht_p=0.8, ADO_p = 0.2, FP_p = 3.2e-5):
    """ Main function to simulate single cell reference data set """

    """ Input:
          source (string): input sequence used to simulate single cells.
          n_snv (int): total number of SNVs to simulate per cell.
          n_sc (int): number of single cells.
          snv_loc (int list): optional list of SNV locations to simulate.
          out (string): output file name, default is 'output'.
          ADO_p (float): allelic dropout rate, must be in [0,1].
          FP_p (float): false positive rate, must be in [0,1].
          ...
          SNV simulation parameters: 
                ts_tv_p (float): transition to transversion rate, must be in [0,1].
                ts_hm_ht_p (float): homozygous to heterozygous transition rate, must be in [0,1].
                AGCT_AGTC_p (float): AG-CT to AG-TC transversion rate, must be in [0,1].
                tv_hm_ht_p (float): homozygous to heterozygous transversion rate, must be in [0,1].
        Output:
          The functions generates fasta files, each containing the simulated
          reference sequence for a given single cell allele.
    """

    # source sequence length
    N = len(source)

    if snv_loc is None:
        # sample snv locations across source
        # NB: avoid extremities of the sequence
        offset = int(0.1*N)
        snv_loc = np.linspace(start=offset, stop=N-offset, num=n_snv).astype(int)
    else:
        snv_loc = snv_loc

    # simulate SNV flags matrix
    (snv_mtx, p_lst) = sim_snv_flags_mtx(n_snv, n_sc)

    # generate alternate genotypes for all snv locations
    T_snv = gen_snv_transition_matrix()
    alt_gtype = {loc: sim_alt_gtype(ref=2*source[loc], T_snv=T_snv)
                 for loc in snv_loc}

    # loop through single cells
    start_time = time.clock() # start timer
    sc_fnames = []
    for c in range(n_sc):
        logger.info('Simulating single cell %d ...'%(c+1))

        sc_gtypes = simulate_snv(source, snv_loc[snv_mtx[:,c]], alt_gtype)

        logger.info('SNV simulation done! Writing sequence to fasta file ...')

        # write simulated single cell sequence to fasta
        outname = out+'_sc%d'%(c+1)
        sc_fnames.append(outname) # store file names
        write_sc(sc_obj = sc_gtypes, seq_type = 'alt_gtype', out = 'before_'+outname) # write sequence before WGA
        write_sc(sc_obj = sc_gtypes, seq_type = 'wga_gtype', out = 'after_'+outname) # write sequence after WGA

        # export SNV locations to bed files
        bed_df = pd.DataFrame({'chrom': "ref_source",
                               'chromStart': snv_loc[snv_mtx[:,c]],
                               'chromEnd': snv_loc[snv_mtx[:,c]]+1})
        bed_df.to_csv(path_or_buf=os.path.join(RES_DIR,'snv_sc%d.bed'%(c+1)), sep='\t', header=False, mode='w', index=False)

        # export ADO locations to bed files

        # export FP locations to bed files


    sim_time = time.clock() - start_time # store simulation time for display
    logger.info('All done! Simulated %d SNVs in %d single cells in %.2f seconds.'%(n_snv,n_sc,sim_time))

    # pack experiment results
    exp_res = dict(sc_fnames=sc_fnames, snv_mtx = snv_mtx, snv_loc=snv_loc)
    return(exp_res)

if __name__ == '__main__':

    # set simulation default parameters
    n_sc_lst = [10, 15, 20]   # number of single cells to simulate
    n_snv = 1000   # True SNV count
    tg_size = 10**6     # size (in bp) of target region 
    start_region = 10**5      # starting sequence index for filtering
    source_name = "/data/chr20.fa"  # source file name
    out_name = "output"
    seed = 10696 # seed for simulation reproducibility
    _ADO_p = 0.2 # allelic dropout rate
    _FP_p = 3.2e-5 # false positive rate

    ## start processing user input ##
    argc = len(sys.argv)
    i = 1
    while(i < argc):
        if(sys.argv[i] == '-e'):
            n_exp = int(sys.argv[i+1]) # number of experiments
        elif(sys.argv[i] == '-n'):
            # number of single cells per experiment
            n_sc_lst = []
            [n_sc_lst.append(int(sys.argv[j])) for j in range(i+1, i+1+n_exp)] 
            i += n_exp - 1
        elif(sys.argv[i] == '-s'):
            # number of SNVs to simulate
            n_snv = int(sys.argv[i+1])
        elif(sys.argv[i] == '-f'):
            # source (input) file name
            source_name = sys.argv[i+1]
        elif(sys.argv[i] == '-ts'):
            # size (in bp) of target region to consider
            tg_size = int(sys.argv[i+1])
        elif(sys.argv[i] == '-r'):
            # starting sequence region to consider 
            start_region = int(sys.argv[i+1])
        elif(sys.argv[i] == '-sd'):
            # seed for simulation reproducibility
            seed = int(sys.argv[i+1])
        elif(sys.argv[i] == '-o'):
            # output file name
            out_name = sys.argv[i+1]
        elif(sys.argv[i] == '-d'):
            # set results directory
            RES_DIR = sys.argv[i+1]
            logger.debug('Writing results in %s...'%RES_DIR)
        elif(sys.argv[i] == '-ado'):
            # set allelic dropout rate
            _ADO_p = sys.argv[i+1]
        elif(sys.argv[i] == '-fp'):
            # set false positive rate
            _FP_p = sys.argv[i+1]
           
        i += 2         # increment counter
    ## end of input processing ##
        
    np.random.seed(seed) # set environment random seed
    
    # read and filter source to target size
    source = filter_source(source=source_name, target_size = tg_size)
    # write filtered source to FASTA
    write_fasta(seq=source, desc='ref_source')
    
    fnames = [] # initialize file names list
    # loop through experiments
    for e in range(len(n_sc_lst)):
        # get number of cells to simulate in current experiment
        n_sc = n_sc_lst[e]
        logger.info("Starting Experiment %d: simulating %d single cell sequences..."%((e+1),n_sc))

        # simulate single cells and write to fasta files (one per allele)
        exp_res = simulate_sc(source=source, n_snv=n_snv, n_sc=n_sc, out=out_name+'_exp%d'%(e+1), ADO_p = _ADO_p, FP_p = _FP_p)

        # store filenames
        fnames.extend(exp_res['sc_fnames'])
        # write SNV flags matrix to file
        write_snv_mtx(snv_mtx = exp_res['snv_mtx'], snv_loc=exp_res['snv_loc'], e=e)


    # write filenames to text file
    ffnames = open(os.path.join(RES_DIR,'ref_fnames.txt'), "w")
    [ffnames.write("after_"+str(i)+"_a1.fa\n"+ "after_"+str(i)+"_a2.fa\n") for i in fnames]
    [ffnames.write("before_"+str(i)+"_a1.fa\n"+"before_"+str(i)+"_a2.fa\n") for i in fnames]
    ffnames.write("ref_source.fa\n")
    ffnames.close()
        
