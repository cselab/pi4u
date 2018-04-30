CC = gcc
CPP = g++ 

CFLAGS = -O3 -std=c99 -Wall 

TARGETS = display_gen display_prop

all: $(TARGETS)

display_gen: display_gen.c gnuplot_i.c
	$(CC) $(CFLAGS) -o display_gen display_gen.c gnuplot_i.c $(LDFLAGS)

display_prop: display_prop.c gnuplot_i.c
	$(CC) $(CFLAGS) -o display_prop display_prop.c gnuplot_i.c $(LDFLAGS)

clean:
	rm -f display_gen display_prop *.o
