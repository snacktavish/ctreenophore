CC = gcc
include_dir = /home/ejmctavish/projects/libpll-1.0.0/src
CFLAGS = -Wall -O2 -fomit-frame-pointer
ARCH = sse3

all: likelihood

likelihood: likelihood.c
	$(CC) $+ -o $@ $(CFLAGS) -I $(include_dir) -L $(include_dir) -lpll-$(ARCH) -lm

clean:
	rm likelihood
