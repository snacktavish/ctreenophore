all: example

example: example.c
	gcc example.c -lpll-sse3 -lm -o example

toy: toy.c
	gcc -g toy.c -lpll-sse3 -lm -o toy

newick_util: newick_util.c
	gcc -g newick_util.c -lpll-sse3 -lm -I/home/ejmctavish/projects/libpll-1.0.0/src -o newick_util


clean: 
	  rm -rf *.o
