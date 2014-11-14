# Letting the user have PLL_HOME in a non-standard
#    location is nice because some folks might want multiple installations 
#    of PLL. I use this sort of thing a lot to have a fully optimization 
#    build and a debugging build of a library.
# The ?= syntax specifies the default if the var is not in the environment.
PLL_HOME ?= /usr/local

# CPPFLAGS are arguments to the C preprocessor - typically -D for defines
#    and -I for includes.
# The := appends, which is nice if you want to be able to tweak CFLAGS
#    in your environment without touching the Makefile
CPPFLAGS := -I${PLL_HOME}/include

# CFLAGS are flags to the C compiler. The distinction between the preprocessor
#  and the compiler often does not matter. But it is possible to just 
#  preprocess with gcc -E
#  my rule of thumb is that if it won't affect the results of that
#  invocation then it belongs in CFLAGS not CPPFLAGS.
#  -g adds "symbols" to your executable, which makes it easier to debug.
CFLAGS := -g

# LDFLAGS are flags to the linker... -L (capital matters!) is the 
#   set of (non-standard) directories to search for libraries in.
LDFLAGS = -L${PLL_HOME}/lib

# LDLIBS is the list of libraries that your executable needs in order to link.
LDLIBS = -lpll-sse3 -lm

all: example

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
example: example.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o example

toy: toy.o
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o toy

newick_util: newick_util.c
	gcc -g newick_util.c -lpll-sse3 -lm -I/home/ejmctavish/projects/libpll-1.0.0/src -o newick_util


clean: 
	  rm -rf *.o
