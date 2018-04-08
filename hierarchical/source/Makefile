#==============================================================================
#	PATHS AND SOURCE FILES
#==============================================================================
GSLROOT=/usr/local

TARGET1 := sample_theta
TARGET2 := sample_theta_fast
TARGET3 := sample_psi
TARGET4 := sample_posterior_theta_fast

COMMON_SRC =  engine_tmcmc.c  tmcmc_stats.c  tmcmc_db.c  auxil.c  spawner.c

SRC1 := $(COMMON_SRC)	likelihoods/fitfun_theta.c		likelihoods/loglike_theta.c
SRC2 := $(COMMON_SRC)	likelihoods/fitfun_theta_fast.c	likelihoods/loglike_theta_fast.c
SRC3 := $(COMMON_SRC)	likelihoods/fitfun_psi.c \
						likelihoods/loglike_psi.c \
						likelihoods/priors.c
SRC4 := $(COMMON_SRC)	likelihoods/fitfun_posterior_theta.c \
						likelihoods/loglike_posterior_theta.c \
						likelihoods/loglike_theta_fast.c \
						likelihoods/priors.c

ALL_SRC := $(SRC1) $(SRC2) $(SRC3) $(SRC4)
TARGETS := $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

OBJDIR = .o
OBJ1:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC1)) )
OBJ2:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC2)) )
OBJ3:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC3)) )
OBJ4:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(SRC4)) )

#==============================================================================
#	COMPILERS AND FLAGS
#==============================================================================


CC := mpicc
LD := mpicc

CFLAGS += -O3 -std=c99
CFLAGS += -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += -Wall -Wno-unused-function
CFLAGS += -I${GSLROOT}/include
CFLAGS += `torc_cflags`
#CFLAGS += -I/cluster/home/garampat/usr/include -D_REENTRANT

LDFLAGS += -L${GSLROOT}/lib -lgsl -lgslcblas
LDFLAGS +=  -lm
#LDFLAGS += `torc_libs`
LDFLAGS += -L/cluster/home/garampat/usr/lib -ltorc -lpthread


COMPILE.c = $(CC) $(CFLAGS) -c -o $@
#LINK.o    = $(LD) $(LDFLAGS) $(LDLIBS) -o $@
LINK.o    = $(LD)  -o $@  


$(shell mkdir -p $(OBJDIR) )

#==============================================================================
#	RULES
#==============================================================================
all: $(TARGETS)




sample_theta: $(OBJ1)
	$(LINK.o) $^ $(LDFLAGS)

sample_theta_fast: $(OBJ2)
	$(LINK.o) $^ $(LDFLAGS)

sample_psi: $(OBJ3)
	$(LINK.o) $^ $(LDFLAGS)

sample_posterior_theta_fast: $(OBJ4)
	$(LINK.o) $^ $(LDFLAGS)





$(OBJDIR)/%.o: %.c
	@mkdir -p "$(@D)"
	$(COMPILE.c) $<


ctags:
	@rm .tags .tags1 tags
	@ctags $(ALL_SRC)


print:
	@echo $(OBJ1)


clean:
	rm -rf  $(TARGETS) $(OBJDIR)


clear:
	rm -rf .tags* tags
