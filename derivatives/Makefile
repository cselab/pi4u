CC = mpicc
CPP = mpic++ 
MPIF90 = mpif90

CFLAGS = -O3 -std=c99 -Wall 

CFLAGS += `torc_cflags`
LDFLAGS += `torc_libs` -lpthread -lm 

CFLAGS += -I../pndl/src
LDFLAGS += -lgsl -lgslcblas

TARGETS = sa_deriv fd_deriv fd_grad

all: $(TARGETS)

sa_deriv: sa_deriv.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o sa_deriv sa_deriv.c auxil.o $(LDFLAGS)

fd_deriv: fd_deriv.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o fd_deriv fd_deriv.c auxil.o ../pndl/src/libpndl.a $(LDFLAGS)

fd_grad: fd_grad.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o fd_grad fd_grad.c auxil.o ../pndl/src/libpndl.a $(LDFLAGS)

clean:
	rm -f sa_deriv fd_deriv fd_grad *.o sHbarbar.txt
