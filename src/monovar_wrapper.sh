#!/bin/bash
source activate monovar
#export PATH=$PATH:../../src/monovar/external/samtools
#export PATH=$PATH:../../src/monovar/src
samtools mpileup -BQ0 -d10000 -f $1 -q 40 -b $2 | monovar -p 0.002 -a 0.2 -t 0.05 -m 2 -f $1 -b $2 -o $3
source deactivate
