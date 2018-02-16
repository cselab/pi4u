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
