all: example

example: example.c
	gcc example.c -lpll-sse3 -lm -o example

toy: toy.c
	gcc -g toy.c -lpll-sse3 -lm -o toy

clean: 
	  rm -rf *.o
