CC = gcc
include_dir = /home/localuser/Library/phylogenetic-likelihood-library/src
CFLAGS = -Wall -O2 -fomit-frame-pointer
ARCH = SSE3

all: likelihood

likelihood: likelihood.c
	$(CC) $+ -o $@ $(CFLAGS) -I $(include_dir) -L $(include_dir) -lpll-$(ARCH) -lm

clean:
	rm likelihood
