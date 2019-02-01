

.PHONY: all build sh

UID = $(shell id -u)
GID = $(shell id -g)

all : build
	cd example && make all

build:
	docker build -t scsim .

env:
	conda env create -f envs/scsim-env.yml

sh:
	docker run --rm -it -v $(PWD):/mnt -w /mnt --user $(UID):$(GID) scsim "source activate scsim"

# vim:ft=make
#
