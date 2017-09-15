# Python

PYTHON_CFLAGS = `python2.7-config --cflags`
PYTHON_LDFLAGS = `python2.7-config --ldflags`
PYTHON_CFLAGS += -I/usr/local/lib/python2.7/site-packages/numpy/core/include

#PYTHON_CFLAGS = `python3.6-config --cflags`
#PYTHON_LDFLAGS = `python3.6-config --ldflags`
#PYTHON_CFLAGS += -I/usr/local/lib/python3.6/site-packages/numpy/core/include


# R
R_BASE=/usr/local/lib/R
#R_BASE=/usr/lib64/R
R_CFLAGS += -g -I$(R_BASE)/include -std=c99
R_LDFLAGS += -L$(R_BASE)/lib -lR #-lRbla

# Compiler

CC = gcc # mpicc
CPP = g++ # mpic++ 
MPIF90 = mpif90

CFLAGS += -O0 -g -ggdb
#CFLAGS = -O3
CFLAGS += -Wall -fPIC


all:Python #Matlab Rlib

Python:libtmcmc.so 
Matlab:mex_tmcmc.mex64
Rlib: libtmcmcR.so

all: $(TARGETS)


UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	mex_suffix = mexa64
endif
ifeq ($(UNAME_S),Darwin)
	mex_suffix = mexmaci64
endif

#GSL_OBJS=../gsl/gsl-2.4/*/.libs/*.o
GSL_OBJS=-lgsl -lgslcblas

libtmcmc.so: engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c fitfun_python.c
	$(CC) $(CFLAGS) $(PYTHON_CFLAGS) -shared -o libtmcmc.so engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c fitfun_python.c $(PYTHON_LDFLAGS) $(GSL_OBJS)
	cp -f libtmcmc.so ../lib_python

mex_tmcmc.mex64: engine_tmcmc.c engine_tmcmc.h tmcmc_stats.c tmcmc_db.c auxil.c fitfun_matlab.c
	mex -O CFLAGS="\$(CFLAGS) -std=c99 -I/usr/local/include" -output mex_tmcmc.$(mex_suffix) engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c fitfun_matlab.c $(GSL_OBJS)
	cp -f mex_tmcmc.$(mex_suffix) ../lib_matlab

libtmcmcR.so: engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c fitfun_R.c
	$(CC) $(CFLAGS) $(R_CFLAGS) -shared -o libtmcmcR.so engine_tmcmc.c tmcmc_stats.c tmcmc_db.c auxil.c fitfun_R.c $(R_LDFLAGS) $(GSL_OBJS)
	cp -f libtmcmcR.so ../lib_R

engine_tmcmc.o: engine_tmcmc.c engine_tmcmc.h
	$(CC) $(CFLAGS) $(PYTHON_CFLAGS) -c engine_tmcmc.c

tmcmc_db.o: tmcmc_db.c
	$(CC) $(CFLAGS) $(PYTHON_CFLAGS) -c tmcmc_db.c

auxil.o: auxil.c
	$(CC) $(CFLAGS) $(PYTHON_CFLAGS) -c auxil.c

fitfun_python.o: fitfun_python.c
	$(CC) $(CFLAGS) $(PYTHON_CFLAGS) -c fitfun_python.c

#python.o: python.c
#	$(CC) $(CFLAGS) $(PYTHON_CFLAGS) -c python.c

clean:
	rm -f libtmcmc.so *.o *.pyc 
	rm -rf libtmcmc.so.dSYM
	rm -f mex_tmcmc.mexa64 mex_tmcmc.mexmaci64
	rm -f libtmcmcR.so
	rm -rf libtmcmcR.so.dSYM
