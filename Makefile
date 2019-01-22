#
# Makefile
# flaherty, 2019-01-22 13:19
#

.PHONY: all build sh

all:
	@echo "Makefile needs your attention"

build :
	docker build -t scsim .

sh :
	docker run --rm -it -v $(PWD):/mnt -w /mnt scsim

# vim:ft=make
#
