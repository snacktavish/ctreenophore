CC = gcc
ARCH = sse3
all: likelihood mini

clean:
	rm -rf *.o likelihood mini
# Letting the user have PLL_HOME in a non-standard
#    location is nice because some folks might want multiple installations 
#    of PLL. I use this sort of thing a lot to have a fully optimization 
#    build and a debugging build of a library.
# The ?= syntax specifies the default if the var is not in the environment.
PLL_HOME ?= /usr/local
include_dir = /home/ejmctavish/projects/libpll-1.0.0/src

# CPPFLAGS are arguments to the C preprocessor - typically -D for defines
#    and -I for includes.
# The := appends, which is nice if you want to be able to tweak CFLAGS
#    in your environment without touching the Makefile
CPPFLAGS := -I${PLL_HOME}/include/pll -I${PLL_HOME}/include -I$(include_dir)

#  -g adds "symbols" to your executable, which makes it easier to debug.
CFLAGS := -Wall -O0 -fomit-frame-pointer -g

# LDFLAGS are flags to the linker... -L (capital matters!) is the 
#   set of (non-standard) directories to search for libraries in.
LDFLAGS = -L${PLL_HOME}/lib -L$(include_dir)

# LDLIBS is the list of libraries that your executable needs in order to link.
LDLIBS = -lpll-$(ARCH) -lm

# This is the implicit rule for making a .o (object) file from .c
# you don't have to include it, it just makes it clear
# which variables (CC, CPPFLAGS, and CFLAGS) are used
# the $@ $< weirdness is using pattern substitution to
# get the name of the target (something.o) associated with the -o flag
# and getting the prerequisite of the rule (something.c) as an arg to
# the compiler.
%.o:%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

# expressing the rule to make the executable in terms of the .o files (like
#    I've done here) is silly in this tiny example. But in general it is 
#    nice because it lets you specify rules for the compilation of each
#    piece into .o precisely, and then describe how to link those .o files
#    into an executable
likelihood: likelihood.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o likelihood

mini: mini.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o mini


