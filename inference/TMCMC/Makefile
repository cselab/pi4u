CC = mpicc
CPP = mpic++ 
MPIF90 = mpif90

CFLAGS = -O3 -std=c99
CFLAGS += -Wall 

LDFLAGS += -lgsl -lgslcblas

CFLAGS += `torc_cflags`
LDFLAGS += `torc_libs` -lpthread -lm 

TARGETS = engine_tmcmc 

all: $(TARGETS)

engine_tmcmc: engine_tmcmc.c engine_tmcmc.h tmcmc_stats.c tmcmc_db.c auxil.c fitfun.c
	$(CC) $(CFLAGS) -o engine_tmcmc engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c $(LDFLAGS)

clean:
	rm -f engine_tmcmc  *.o

clear:
	rm -rf curgen_db*.txt curres_db*.txt runinfo.txt full_db*.txt fitness.txt final.txt
