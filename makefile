CC=gcc
CFLAGS=-Wall -g -std=c99

all: hydro-loops

hydro-loops: hydro-loops_functions.o hydro-loops_main.o  
	$(CC) $(CFLAGS) hydro-loops_functions.o hydro-loops_main.o -o hydro-loops -lm

hydro-loops_functions.o: hydro-loops_functions.c
	$(CC) -c  $(CFLAGS) hydro-loops_functions.c 

hydro-loops_main.o: hydro-loops_main.c
	$(CC) -c  $(CFLAGS) hydro-loops_main.c

clean: 
	rm hydro-loops*.o hydro-loops
