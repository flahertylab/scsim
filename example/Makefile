#
# Makefile
# flaherty, 2019-01-22 15:55
#

.PHONY: all clean
SNAKEMAKE=docker run --rm -v $(PWD):/mnt -w /mnt scsim snakemake

all:
	$(SNAKEMAKE) all

clean :
	$(SNAKEMAKE) clean

# vim:ft=make
#
