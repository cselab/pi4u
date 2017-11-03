CC = mpicc
CPP = mpic++
MPIF90 = mpif90


CFLAGS = -O3 -std=c99
CFLAGS += -Wall -D_XOPEN_SOURCE=700 -D_BSD_SOURCE


GSLROOT=/usr/local
CFLAGS  += -I${GSLROOT}/include


LDFLAGS += -L${GSLROOT}/lib -lgsl -lgslcblas

CFLAGS += `torc_cflags`
LDFLAGS += `torc_libs` -lpthread -lm

TARGETS = sample_theta  sample_theta_fast  sample_psi

all: $(TARGETS)


sample_theta: engine_tmcmc.c engine_tmcmc.h spawner.h tmcmc_stats.c tmcmc_db.c auxil.c fitfun.h
	$(CC) $(CFLAGS) -o sample_theta engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c spawner.c fitfun_theta.c $(LDFLAGS)


sample_theta_fast: engine_tmcmc.c engine_tmcmc.h spawner.h tmcmc_stats.c tmcmc_db.c auxil.c fitfun.h
	$(CC) $(CFLAGS) -o sample_theta_fast engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c spawner.c fitfun_theta_fast.c $(LDFLAGS)


sample_psi: engine_tmcmc.c engine_tmcmc.h spawner.h tmcmc_stats.c tmcmc_db.c auxil.c fitfun.h
	$(CC) $(CFLAGS) -o sample_psi   engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c spawner.c fitfun_psi.c $(LDFLAGS)


clean:
	rm -f sample_theta sample_theta_fast sample_psi  *.o

clear:
	rm -rf *.txt
