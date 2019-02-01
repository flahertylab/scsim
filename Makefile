

.PHONY: all build sh

UID = $(shell id -u)
GID = $(shell id -g)

all : build
	cd example && make all

build:
	docker build -t scsim .

sh:
	docker run --rm -it -v $(PWD):/mnt -w /mnt --user $(UID):$(GID) scsim

# vim:ft=make
#
