#==============================================================================
#	PATHS AND SOURCE FILES
#==============================================================================
GSLROOT=$(HOME)/usr/gsl

TARGET1 := sample_theta
TARGET2 := sample_theta_fast
TARGET3 := sample_psi
TARGET4 := sample_posterior_theta
TARGET5 := optimize_theta
TARGET6 := optimize_theta_fast
TARGET7 := dram_theta_fast
TARGET8 := dram_psi


SRC := ../source

TMCMC_SRC =  	$(SRC)/TMCMC/tmcmc_run.c \
				$(SRC)/TMCMC/tmcmc_engine.c \
				$(SRC)/TMCMC/tmcmc_stats.c \
				$(SRC)/TMCMC/tmcmc_db.c \
				$(SRC)/TMCMC/tmcmc_aux.c \


CMAES_SRC = $(SRC)/CMAES/engine_cmaes.c \
			$(SRC)/CMAES/cmaes.c

DRAM_SRC = 	$(SRC)/DRAM/run_dram.c \
			$(SRC)/DRAM/dram.c
						
PRIORS = 	$(SRC)/priors/priors.c \
			$(SRC)/priors/myrand.c


SRC1 := $(TMCMC_SRC) $(PRIORS)		$(SRC)/likelihoods/fitfun_theta.c \
									$(SRC)/likelihoods/loglike_theta.c \
									$(SRC)/likelihoods/spawner.c

SRC2 := $(TMCMC_SRC) $(PRIORS)		$(SRC)/likelihoods/fitfun_theta_fast.c \
									$(SRC)/likelihoods/loglike_theta_fast.c

SRC3 := $(TMCMC_SRC) $(PRIORS)		$(SRC)/likelihoods/fitfun_psi.c \
									$(SRC)/likelihoods/loglike_psi.c \
									$(SRC)/likelihoods/priors.c

SRC4 := $(TMCMC_SRC) $(PRIORS)		$(SRC)/likelihoods/fitfun_posterior_theta.c \
									$(SRC)/likelihoods/loglike_posterior_theta.c \
									$(SRC)/likelihoods/loglike_theta_fast.c \
									$(SRC)/likelihoods/priors.c

SRC5 := $(CMAES_SRC) $(PRIORS)    	$(SRC)/likelihoods/fitfun_theta.c \
									$(SRC)/likelihoods/loglike_theta.c \
									$(SRC)/likelihoods/spawner.c

SRC6 := $(CMAES_SRC) $(PRIORS)		$(SRC)/likelihoods/fitfun_theta_fast.c \
									$(SRC)/likelihoods/loglike_theta_fast.c\

SRC7 := $(DRAM_SRC) $(PRIORS)		$(SRC)/likelihoods/fitfun_theta_fast.c \
									$(SRC)/likelihoods/loglike_theta_fast.c 

SRC8 := $(DRAM_SRC)	$(PRIORS)		$(SRC)/likelihoods/fitfun_psi.c \
									$(SRC)/likelihoods/loglike_psi.c \
									$(SRC)/likelihoods/priors.c

ALL_SRC := $(SRC1) $(SRC2) $(SRC3) $(SRC4) $(SRC5)  $(SRC6)  $(SRC7)  $(SRC8)
TARGETS := $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7) $(TARGET8)

OBJDIR = .o
TMP=$(subst ../,,$(SRC1))
OBJ1:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)) )
TMP=$(subst ../,,$(SRC2))
OBJ2:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)) )
TMP=$(subst ../,,$(SRC3))
OBJ3:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)) )
TMP=$(subst ../,,$(SRC4))
OBJ4:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)) )
TMP=$(subst ../,,$(SRC5))
OBJ5:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)) )
TMP=$(subst ../,,$(SRC6))
OBJ6:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)) )
TMP=$(subst ../,,$(SRC7))
OBJ7:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)) )
TMP=$(subst ../,,$(SRC8))
OBJ8:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)) )

#==============================================================================
#	COMPILERS AND FLAGS
#==============================================================================
use_torc?=0
use_omp?=1

CC := gcc
LD := gcc

ifeq ($(use_torc),1)
	CC := mpicc
	LD := mpicc
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
CFLAGS += -I./$(SRC)/TMCMC

LDLIBS += -L${GSLROOT}/lib -lgsl -lgslcblas -lm -lpthread


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

dram_psi: $(OBJ8)
	$(LINK.o) $^ $(LDLIBS)


$(OBJDIR)/%.o: ../%.c
	@mkdir -p "$(@D)"
	$(COMPILE.c) $<


ctags:
	@rm .tags .tags1 tags
	@ctags $(ALL_SRC)


print:
	@echo $(OBJ6)


clean:
	rm -rf  $(TARGETS) $(OBJDIR)


clear:
	rm -rf .tags* tags
