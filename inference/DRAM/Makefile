CC=gcc
CFLAGS=-O3 -Wall -std=c99 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
#CFLAGS+= -DDEBUG=1
LDFLAGS=-lgsl -lgslcblas -lm

all: normaltest bananatest generaltest

normaltest: normaltest.c dramrun.c
	$(CC) $(CFLAGS) -o normaltest normaltest.c $(LDFLAGS)

bananatest: bananatest.c dramrun.c
	$(CC) $(CFLAGS) -o bananatest bananatest.c $(LDFLAGS)

generaltest: generaltest.c dramrun.c fitfun.c
	$(CC) $(CFLAGS) -o generaltest generaltest.c fitfun.c $(LDFLAGS)

clean:
	rm -f normaltest bananatest generaltest

clear:
	rm -f *test_chain.txt
