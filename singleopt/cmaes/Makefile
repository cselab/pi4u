use_torc?=0

CC=gcc
CFLAGS=-O3 -std=c99 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
LDFLAGS=-lm -lpthread

ifeq ($(use_torc),1)
CC = mpicc
CFLAGS += -D_USE_TORC_=1 `torc_cflags`
LDFLAGS += `torc_libs`
endif

ifeq ($(use_torc),1)
TARGETS = engine_cmaes engine_cmaes_async
else
TARGETS = engine_cmaes engine_cmaes_async_note
endif

all: $(TARGETS)

engine_cmaes: engine_cmaes.c cmaes.c fitfun.c
	$(CC) $(CFLAGS) -o engine_cmaes engine_cmaes.c cmaes.c $(LDFLAGS)

engine_cmaes_async: engine_cmaes.c cmaes.c fitfun.c
	$(CC) $(CFLAGS) -o engine_cmaes_async engine_cmaes_async.c cmaes.c $(LDFLAGS)

.PHONY: engine_cmaes_async_note
engine_cmaes_async_note:
	@echo "To build the experimental async version of CMAES you must use: make use_torc=1"

.PHONY: clean
clean:
	rm -f engine_cmaes engine_cmaes_async *.o

.PHONY: clear
clear:
	rm -rf *.dat curgen_db*.txt actparcmaes.par
