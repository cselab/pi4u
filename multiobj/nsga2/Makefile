use_torc?=0
use_omp?=0

# Makefile for compiling NSGA-II source code
CC=gcc
LD=gcc
#CFLAGS=-Wall -ansi -pedantic -g

CFLAGS = -g -O3 -std=c99 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += -Wall
CFLAGS += -DIODUMP=1 -DBINARY=0 
LDFLAGS += -pthread -lm

ifeq ($(use_torc),1)
CC = mpicc
LD = mpicc
CFLAGS += -D_USE_TORC_=1 `torc_cflags` #-I/users/sverma/usr/torc/include
LIBS += `torc_libs` #-L/users/sverma/usr/torc/lib -ltorc
endif

ifeq ($(use_omp),1)
CFLAGS += -D_USE_OPENMP_=1 -fopenmp
LIBS += -fopenmp
endif

RM=rm -f

OBJS1:=allocate.o auxiliary.o crossover.o crowddist.o decode.o display.o dominance.o eval.o fillnds.o fitfun_inline.o initialize.o list.o merge.o mutation.o nsga2r.o problemdef.o rand.o rank.o report.o sort.o tourselect.o
MAIN1=nsga2r

OBJS2:=allocate.o auxiliary.o crossover.o crowddist.o decode.o display.o dominance.o eval.o fillnds.o fitfun_external.o initialize.o list.o merge.o mutation.o nsga2r.o problemdef.o rand.o rank.o report.o sort.o tourselect.o
MAIN2=nsga2r_ext

all:$(MAIN1) $(MAIN2)

$(MAIN1):$(OBJS1)
	$(LD) $(LDFLAGS) $(OBJS1) -o $(MAIN1) $(LIBS)

$(MAIN2):$(OBJS2)
	$(LD) $(LDFLAGS) $(OBJS2) -o $(MAIN2) $(LIBS)

%.o: %.c global.h rand.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) $(OBJS1) $(MAIN1) $(OBJS2) $(MAIN2)

clear:
	$(RM) eval_db_*.txt curgen_db_*.txt *.out
