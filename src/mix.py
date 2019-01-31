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

def mix_fastq(jsonconfig):
      # Load experiment configuration from json file

      mix_config = json.loads(jsonconfig)

      # Set the rng
      np.random.seed(mix_config['rng_seed'])

      # Remove output files if they exist
      remove_file(mix_config['mix_path'])

      # read-sample-write each single cell
      for i in range(1, len(mix_config) - 2):
            logger.info("Processing single cell {}...".format(i))
            cell_readcount = int(np.floor(mix_config['total_reads'] * mix_config['sc{}'.format(i)]['fraction']))
            etl_mix(source_path = mix_config['sc{}'.format(i)]['path'],
                  readcount = cell_readcount,
                  mix_path = mix_config['mix_path'])


def etl_mix(source_path, readcount, mix_path):
    """Load reads from source fq.gz file, sample readcount reads, and append to mix fq."""

    logger.info("Loading reads...")
    reads = load_reads(source_path)

    # Check that readcount is less than total number of reads available
    if(readcount > len(reads)):
        sys.exit(1)

    # Sample reads
    logger.info("Sampling reads...")
    sample_idx = np.random.choice(len(reads), size=readcount, replace=False)

    # Write reads
    write_reads(reads, sample_idx, mix_path)

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

