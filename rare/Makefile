CPP = mpic++ 

CFLAGS = -O3 -Wall

CFLAGS += `torc_cflags`
LDFLAGS += `torc_libs` -lpthread -lm 

# GSL
LDFLAGS += -lgsl -lgslcblas

# PRIMARY

TARGETS = engine_ss1 engine_ss2 

all: $(TARGETS)

engine_ss1: engine_ss1.cpp ss_db.cpp auxil.o fitfun.c
	$(CPP) $(CFLAGS) -o engine_ss1 engine_ss1.cpp ss_db.cpp auxil.o $(LDFLAGS)

engine_ss2: engine_ss2.cpp ss_db.cpp auxil.o fitfun.c
	$(CPP) $(CFLAGS) -o engine_ss2 engine_ss2.cpp ss_db.cpp auxil.o $(LDFLAGS)

clean:
	rm -f engine_ss1 engine_ss2 *.o

clear:
	rm -rf samples*.txt seeds*.txt Pc.txt

