all: example

example: example.c
	gcc example.c -lpll-sse3 -lm -o example
clean: 
	  rm -rf *.o
