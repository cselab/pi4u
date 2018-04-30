#==============================================================================
#	PATHS AND SOURCE FILES
#==============================================================================
GSLROOT=/usr/local

TARGET1 := engine_tmcmc


TMCMC := .


TMCMC_SRC =  	$(TMCMC)/tmcmc_run.c \
				$(TMCMC)/tmcmc_stats.c \
				$(TMCMC)/tmcmc_db.c \
				$(TMCMC)/tmcmc_aux.c \
				$(TMCMC)/tmcmc_engine.c \

PRIORS	= 	../priors/priors.c\
			../priors/myrand.c

SRC1 := $(TMCMC_SRC) $(PRIORS) fitfun.c


ALL_SRC := $(SRC1) 
TARGETS := $(TARGET1)

OBJDIR = .o
TMP=$(subst ../,,$(SRC1))
OBJ1:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)))
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


engine_tmcmc: $(OBJ1)
	$(LINK.o) $^ $(LDLIBS)




$(OBJDIR)/%.o: %.c
	@mkdir -p "$(@D)"
	$(COMPILE.c) $<


$(OBJDIR)/priors/%.o: ../priors/%.c
	@mkdir -p "$(@D)"
	$(COMPILE.c) $<


tmp:
	@echo $(OBJ1)


ctags:
	@rm -f .tags .tags1 tags
	@ctags $(ALL_SRC)  ../priors/*.c  ../priors/*.h


print:


clean:
	rm -rf  $(TARGETS) $(OBJDIR)


clear:
	rm -rf .tags* tags *.txt
