


.PHONY: all clean
UID = $(shell id -u)
GID = $(shell id -g)
SNAKEMAKE=docker run --rm -v $(PWD):/mnt -w /mnt --user $(UID):$(GID) scsim snakemake --use-conda

all:
	$(SNAKEMAKE) all

clean :
	$(SNAKEMAKE) clean

# vim:ft=make
#
