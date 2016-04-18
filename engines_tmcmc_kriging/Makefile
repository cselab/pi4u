CC = cc
CPP = cc 
MPIF90 = ftn

#CFLAGS = -O0 -g -ggdb -I../pndl/src
CFLAGS = -g -ggdb -O3 #-fopenmp #-std=c99
CFLAGS += -Wall #-Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-result -Wall # -w #-Wall

# TORC
#CFLAGS += `torc_cflags`
#LDFLAGS += `torc_libs` -lpthread -lm 
CFLAGS += -I/users/lina/usr/include
LDFLAGS += -L/users/lina/usr/lib -ltorc -lpthread -lm 

# GSL
LDFLAGS += -lgsl -lgslcblas

TARGETS = display_gen engine_tmcmc # engine_cmaes
#TARGETS += engine_tmcmc_sn
#TARGETS += sa_deriv fd_deriv fd_grad
#TARGETS += engine_ss1 engine_ss2 
TARGETS += propagation_tool

# OpenANN
#CFLAGS += -I${HOME}/usr/fann/include
#LDFLAGS += ../surrogates/OpenANN/nnet.o `pkg-config --cflags --libs openann` -lalglib -lstdc++

# FLANN
#CFLAGS += -I/users/lina/usr/flann/include/flann -I/users/lina/usr/flann/include
#LDFLAGS += -L/users/lina/usr/flann/lib -lflann


all: $(TARGETS)

engine_tmcmc: engine_tmcmc.c engine_tmcmc.h tmcmc_stats.c tmcmc_db.c gnuplot_i.c auxil.c fitfun.c tmcmc_kriging.c spawner.h
	$(CC) $(CFLAGS) -o engine_tmcmc engine_tmcmc.c tmcmc_stats.c tmcmc_db.c gnuplot_i.c auxil.c tmcmc_kriging.c spawner.c $(LDFLAGS)

display_gen: display_gen.c gnuplot_i.c
	$(CC) $(CFLAGS) -o display_gen display_gen.c gnuplot_i.c $(LDFLAGS)

engine_cmaes: engine_cmaes.c cmaes.c fitfun.c
	$(CC) $(CFLAGS) -o engine_cmaes engine_cmaes.c cmaes.c auxil.c $(LDFLAGS)

engine_ss1: engine_ss1.cpp ss_db.cpp auxil.o fitfun.c
	$(CPP) $(CFLAGS) -o engine_ss1 engine_ss1.cpp ss_db.cpp auxil.o $(LDFLAGS)

engine_ss2: engine_ss2.cpp ss_db.cpp auxil.o fitfun.c
	$(CPP) $(CFLAGS) -o engine_ss2 engine_ss2.cpp ss_db.cpp auxil.o $(LDFLAGS)

auxil.o: auxil.c
	$(CC) $(CFLAGS) -c auxil.c 

engine_tmcmc_sn: engine_tmcmc_sn.c engine_tmcmc.h tmcmc_stats.c tmcmc_db.c gnuplot_i.c auxil.c fitfun.c
	$(MPIF90) $(CFLAGS) -D_TMCMC_SN_ -o engine_tmcmc_sn engine_tmcmc_sn.c tmcmc_stats.c tmcmc_db.c gnuplot_i.c auxil.c ../pndl/src/libpndl.a $(LDFLAGS)

sa_deriv: sa_deriv.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o sa_deriv sa_deriv.c auxil.o $(LDFLAGS)

fd_deriv: fd_deriv.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o fd_deriv fd_deriv.c auxil.o ../pndl/src/libpndl.a $(LDFLAGS)

fd_grad: fd_grad.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o fd_grad fd_grad.c auxil.o ../pndl/src/libpndl.a $(LDFLAGS)

propagation_tool: propagation_tool.c auxil.c fitfun.c
	$(CC) $(CFLAGS) -o propagation_tool propagation_tool.c auxil.c $(LDFLAGS)

clean:
	rm -f engine_tmcmc display_gen engine_cmaes engine_ss1 engine_ss2 engine_tmcmc_sn sa_deriv fd_deriv fd_grad propagation_tool *.o

clear:
	rm -rf lsf.o* *.dat curgen_db*.txt curres_db*.txt full_db*.txt weights_*.txt P*S*T* runinfo.txt samples*.txt seeds*.txt Pc.txt actparcmaes.par myNN_???.txt scale_myNN_???.txt

backup:
	scp -r $(PWD) chatzidp@bulldozer.ethz.ch:
