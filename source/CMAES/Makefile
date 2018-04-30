GSLROOT:=$(HOME)/usr/gsl/

SRC := 	engine_cmaes.c \
		cmaes.c \
		fitfun.c \
		../priors/priors.c \
		../priors/myrand.c


ALL_SRC := $(SRC)
TARGETS := cma

OBJDIR = .o
TMP=$(subst ../,,$(SRC))
OBJ:= $(patsubst %,$(OBJDIR)/%.o,$(basename  $(TMP)))




#==============================================================================
#	COMPILERS AND FLAGS
#==============================================================================
use_torc?=1

CC := mpicc
LD := mpicc

ifeq ($(use_torc),1)
	CC = mpicc
	CFLAGS += -D_USE_TORC_=1 `torc_cflags`
	LDLIBS += `torc_libs`
endif


CFLAGS += -O3 -std=c99
CFLAGS += -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
CFLAGS += -Wall -Wno-unused-function
CFLAGS += -I${GSLROOT}/include

LDLIBS += -L${GSLROOT}/lib -lgsl -lgslcblas  -lm -lpthread


COMPILE.c = $(CC) $(CFLAGS) -c -o $@
LINK.o    = $(LD) $(LDFLAGS) -o $@


$(shell mkdir -p $(OBJDIR) )

#==============================================================================
#	RULES
#==============================================================================
all: $(TARGETS)


cma: $(OBJ)
	$(LINK.o) $^ $(LDLIBS)


$(OBJDIR)/%.o: %.c
	@mkdir -p "$(@D)"
	$(COMPILE.c) $<


$(OBJDIR)/priors/%.o: ../priors/%.c
	@mkdir -p "$(@D)"
	$(COMPILE.c) $<



ctags:
	@rm .tags .tags1 tags
	@ctags $(ALL_SRC)


print:


clean:
	rm -rf  $(TARGETS) $(OBJDIR)


clear:
	rm -rf .tags* tags *.txt *.dat actparcmaes.par *.err
