CC = mpicc
CPP = mpic++ 
MPIF90 = mpif90

CFLAGS = -O3 -std=c99 -Wall -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += `torc_cflags`

LDFLAGS += `torc_libs` -lpthread -lm 

TARGETS = engine_cmaes engine_cmaes_async

all: $(TARGETS)

engine_cmaes: engine_cmaes.c cmaes.c fitfun.c
	$(CC) $(CFLAGS) -o engine_cmaes engine_cmaes.c cmaes.c $(LDFLAGS)

engine_cmaes_async: engine_cmaes.c cmaes.c fitfun.c
	$(CC) $(CFLAGS) -o engine_cmaes_async engine_cmaes_async.c cmaes.c $(LDFLAGS)

clean:
	rm -f engine_cmaes engine_cmaes_async *.o

clear:
	rm -rf *.dat curgen_db*.txt actparcmaes.par
