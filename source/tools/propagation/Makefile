CC = mpicc
CPP = mpic++ 
MPIF90 = mpif90

CFLAGS = -O3 -std=c99 -Wall 

CFLAGS += `torc_cflags`

# GSL
#LDFLAGS += -lgsl -lgslcblas

# PRIMARY
LDFLAGS += `torc_libs` -lpthread -lm 


TARGETS += propagation_tool sensitivity_tool

all: $(TARGETS)

propagation_tool: propagation_tool.c fitfun.c
	$(CC) $(CFLAGS) -o propagation_tool propagation_tool.c $(LDFLAGS)

sensitivity_tool: sensitivity_tool.c fitfun.c
	$(CC) $(CFLAGS) -o sensitivity_tool sensitivity_tool.c $(LDFLAGS)

clean:
	rm -f propagation_tool sensitivity_tool *.o
