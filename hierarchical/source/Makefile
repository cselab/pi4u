#==============================================================================
#	PATHS AND SOURCE FILES
#==============================================================================
GSLROOT=/usr/local

TARGET1 := sample_theta
TARGET2 := sample_theta_fast
TARGET3 := sample_psi
TARGET4 := sample_posterior_theta
TARGET5 := optimize_theta
TARGET6 := optimize_theta_fast
TARGET7 := dram_theta_fast

TMCMC := TMCMC

TMCMC_SRC =  	$(TMCMC)/engine_tmcmc.c \
				$(TMCMC)/tmcmc_stats.c \
				$(TMCMC)/tmcmc_db.c \
				$(TMCMC)/auxil.c \

CMAES_SRC = CMAES/engine_cmaes.c \
			CMAES/cmaes.c

DRAM_SRC = 	DRAM/dram.c \
			DRAM/auxil.c


SRC1 := $(TMCMC_SRC)	likelihoods/fitfun_theta.c \
						likelihoods/loglike_theta.c \
						likelihoods/spawner.c
SRC2 := $(TMCMC_SRC)	likelihoods/fitfun_theta_fast.c \
						likelihoods/loglike_theta_fast.c
SRC3 := $(TMCMC_SRC)	likelihoods/fitfun_psi.c \
						likelihoods/loglike_psi.c \
						likelihoods/priors.c
SRC4 := $(TMCMC_SRC)	likelihoods/fitfun_posterior_theta.c \
						likelihoods/loglike_posterior_theta.c \
						likelihoods/loglike_theta_fast.c \
						likelihoods/priors.c
SRC5 := $(CMAES_SRC)    likelihoods/fitfun_theta.c \
						likelihoods/loglike_theta.c \
						likelihoods/spawner.c
SRC6 := $(CMAES_SRC)    likelihoods/fitfun_theta_fast.c \
						likelihoods/loglike_theta_fast.c

SRC7 := $(DRAM_SRC)		likelihoods/fitfun_theta_fast.c \
						likelihoods/loglike_theta_fast.c 


ALL_SRC := $(SRC1) $(SRC2) $(SRC3) $(SRC4) $(SRC5)  $(SRC6)  $(SRC7)
TARGETS := $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7)

OBJDIR = .o
OBJ1:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC1)) )
OBJ2:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC2)) )
OBJ3:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC3)) )
OBJ4:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC4)) )
OBJ5:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC5)) )
OBJ6:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC6)) )
OBJ7:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC7)) )

#==============================================================================
#	COMPILERS AND FLAGS
#==============================================================================
use_torc?=1
use_omp?=0

CC := mpicc
LD := mpicc

ifeq ($(use_torc),1)
	CC = mpicc
	CFLAGS += -D_USE_TORC_=1 `torc_cflags`
	LDLIBS += `torc_libs`
endif

ifeq ($(use_omp),1)
	CFLAGS += -D_USE_OPENMP_=1 -fopenmp
	LDLIBS += -fopenmp
endif




CFLAGS += -O3 -std=c99
CFLAGS += -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += -Wall -Wno-unused-function
CFLAGS += -I${GSLROOT}/include
CFLAGS += -I./$(TMCMC)

LDLIBS += -L${GSLROOT}/lib -lgsl -lgslcblas -lm


COMPILE.c = $(CC) $(CFLAGS) -c -o $@
LINK.o    = $(LD) $(LDFLAGS) -o $@


$(shell mkdir -p $(OBJDIR) )

#==============================================================================
#	RULES
#==============================================================================
all: $(TARGETS)


sample_theta: $(OBJ1)
	$(LINK.o) $^ $(LDLIBS)

sample_theta_fast: $(OBJ2)
	$(LINK.o) $^ $(LDLIBS)

sample_psi: $(OBJ3)
	$(LINK.o) $^ $(LDLIBS)

sample_posterior_theta: $(OBJ4)
	$(LINK.o) $^ $(LDLIBS)

sample_posterior_theta_fast: $(OBJ5)
	$(LINK.o) $^ $(LDLIBS)

optimize_theta: $(OBJ5)
	$(LINK.o) $^ $(LDLIBS)

optimize_theta_fast: $(OBJ6)
	$(LINK.o) $^ $(LDLIBS)

dram_theta_fast: $(OBJ7)
	$(LINK.o) $^ $(LDLIBS)


$(OBJDIR)/%.o: %.c
	@mkdir -p "$(@D)"
	$(COMPILE.c) $<


ctags:
	@rm .tags .tags1 tags
	@ctags $(ALL_SRC)


print:


clean:
	rm -rf  $(TARGETS) $(OBJDIR)


clear:
	rm -rf .tags* tags
