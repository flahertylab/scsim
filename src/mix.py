import os
import sys
import resource
import subprocess
import argparse
import logging
import itertools
import collections
import json
import numpy as np
import pdb

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def make_expts(b_fq_toc, b_compbined_alleles, b_r1, b_r2, bulk_toc, G_ij, seed):
    
    BULK_FQ = b_r1 + b_r2  # combined the r1 and r2 lists
    with open(bulk_toc, 'w') as TOC:
        for bulk_sample in range(len(BULK_FQ)//2):
            dist = G_ij[bulk_sample]
            experiment = {'total_reads': 100000, 'rng_seed': seed,
                          'mix_path_r1': BULK_FQ[bulk_sample],
                          'mix_path_r2': BULK_FQ[len(BULK_FQ)//2 + bulk_sample]}
            TOC.write('Bulk sample {} consists of:\n'.format(bulk_sample + 1))
            b_fq_fnames = open(b_fq_toc, 'r').readlines()
            logger.info("mixing fastqs")
            for num in range(len(b_fq_fnames) // 2):
                # num = int(proto.split('/')[-1].split('_')[2][2:])
                proto_r1 = b_fq_fnames[2*num].strip()
                proto_r2 = b_fq_fnames[2*num+1].strip()
                experiment['sc{}'.format(num+1)] = { 'fraction': dist[num], 
                                                    'path_r1': proto_r1,
                                                    'path_r2': proto_r2 }
                TOC.write('    {} of mutant reference {}\n'.format(dist[num - 1], num))
            experiment = json.dumps(experiment)
            mix_fastq(experiment)


def mix_fastq(jsonconfig):
    # Load experiment configuration from json filet
   
    mix_config = json.loads(jsonconfig)
 
    # Set the rng
    np.random.seed(mix_config['rng_seed'])
 
    # Remove output files if they exist
    remove_file(mix_config['mix_path_r1'])
    remove_file(mix_config['mix_path_r2'])
 
    # read-sample-write each single cell
    for i in range(1, len(mix_config) - 3):
          logger.info("Processing single cell {}...".format(i))
          cell_readcount = int(np.floor(mix_config['total_reads'] * mix_config['sc{}'.format(i)]['fraction']))
          etl_mix(source_path_r1 = mix_config['sc{}'.format(i)]['path_r1'],
                  source_path_r2 = mix_config['sc{}'.format(i)]['path_r2'],
                  readcount = cell_readcount,
                  mix_path_r1 = mix_config['mix_path_r1'],
                  mix_path_r2 = mix_config['mix_path_r2'])


def etl_mix(source_path_r1, source_path_r2, readcount, mix_path_r1, mix_path_r2):
    """Load reads from source fq.gz file, sample readcount reads, and append to mix fq."""

    logger.info("Loading reads...")
    read1 = load_reads(source_path_r1)
    read2 = load_reads(source_path_r2)

    # Check that readcount is less than total number of reads available
    if(readcount > len(read1)):
        sys.exit(1)
    if(readcount > len(read2)):
        sys.exit(1)

    # Sample reads
    logger.info("Sampling reads...")
    sample_idx = np.random.choice(len(read1), size=readcount, replace=False)

    # Write reads
    write_reads(read1, sample_idx, mix_path_r1)
    write_reads(read2, sample_idx, mix_path_r2)

def load_reads(source_path):
    """Returns reads from source path fq.gz file."""

    logger.info("Uncompressing and loading {}...".format(source_path))
    readlines = subprocess.check_output(["gunzip", "-c", source_path]).splitlines()

    # Populate a list of lists for the read data
    logger.info("Compiling data structure...")
    read_list = []
    for i in range(0, len(readlines), 4):
        read_id = readlines[i].decode('utf-8')
        read_seq = readlines[i+1].decode('utf-8')
        read_id2 = readlines[i+2].decode('utf-8')
        read_qc = readlines[i+3].decode('utf-8')
        read_list.append([read_id, read_seq, read_id2, read_qc])

    logger.debug("Using {} Mb of memory.".format(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000))

    return read_list

def write_reads(reads, sample_idx, mix_path):

    # Write R1 reads to file
    logger.info("Writing mixture fq...")
    with open(mix_path, 'a') as fh:
        for idx in sample_idx:
            fh.write("\n".join(reads[idx]))
            fh.write("\n")

def remove_file(filename):
      """Deletes a file if it exists."""
      try:
            os.remove(filename)
      except OSError:
            pass

