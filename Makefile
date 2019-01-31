#
# Makefile
# flaherty, 2019-01-22 13:19
#

.PHONY: all build sh

build :
	docker build -t scsim .

sh :
	docker run --rm -it -v $(PWD):/mnt -w /mnt scsim

# vim:ft=make
#
