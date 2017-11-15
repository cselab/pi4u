# Makefile for compiling NSGA-II source code
CC=mpicc
LD=mpicc
RM=rm -f
#CFLAGS=-Wall -ansi -pedantic -g
CFLAGS=-Wall --std=c99 -g -DIODUMP=1 -DBINARY=0 
#CFLAGS+=-fopenmp
CFLAGS+= -D_XOPEN_SOURCE=700 -D_BSD_SOURCE

#LDFLAGS+=-fopenmp
#LIBS+=`torc_libs` -lm

LIBS+= `torc_libs` -lm  # -L/users/sverma/usr/torc/lib -ltorc
CFLAGS+=`torc_cflags` #  -I/users/sverma/usr/torc/include

OBJS:=allocate.o auxiliary.o crossover.o crowddist.o decode.o display.o dominance.o eval.o fillnds.o fitfun.o initialize.o list.o merge.o mutation.o nsga2r.o problemdef.o rand.o rank.o report.o sort.o tourselect.o
MAIN=nsga2r
all:$(MAIN)

$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) $(LIBS)

%.o: %.c global.h rand.h
	$(CC) $(CFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(MAIN)

clear:
	$(RM) eval_db_*.txt curgen_db_*.txt *.out
