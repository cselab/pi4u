CC=gcc
CFLAGS=-O3 -Wall -std=c99 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
LDFLAGS=-lgsl -lgslcblas -lm

all: dram


dram: run_dram.c dram.c fitfun.c auxil.c auxil.h
	$(CC) $(CFLAGS) -o dram run_dram.c dram.c auxil.c fitfun.c $(LDFLAGS)

clean:
	rm -f dram

clear:
	rm -f chain.txt
