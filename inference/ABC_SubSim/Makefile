use_torc?=0

CPP=g++
CPPFLAGS=-std=c++11 -O2 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CC=gcc
CFLAGS=-O2
LDFLAGS=-lgsl -lgslcblas

# uncomment to enable parallel execution
ifeq ($(use_torc),1)
CPP = mpic++
CC = mpicc
CPPFLAGS += -D_USE_TORC_=1 `torc_cflags`
LDFLAGS += `torc_libs`
endif

all: abc_subsim_torc.exe display_gen.exe

abc_subsim_torc.exe: abc_subsim.cpp fitfun.cpp abc_subsim.h problem.h random.h
	${CPP} ${CPPFLAGS} abc_subsim.cpp fitfun.cpp ${LDFLAGS} -o abc_subsim_torc.exe


display_gen.exe: display_gen.c gnuplot_i.c
	$(CC) -o display_gen.exe display_gen.c gnuplot_i.c

clean:
	rm -f *.exe

clear:
	rm -f curgen_db*.txt


cleanall:
	rm -f *.exe
	rm -f curgen_db*.txt
