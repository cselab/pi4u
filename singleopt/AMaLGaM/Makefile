CC = gcc
CFLAGS = #-fopenmp
LIBS = -lm

TCC = mpicc
TCFLAGS = `torc_cflags` 
TLIBS = `torc_libs` -lm

#CFLAGS += -D_RANDOMIZE_
#TCFLAGS += -D_RANDOMIZE_

all: AMaLGaM-Full-Free AMaLGaM-Full-Free_torc AMaLGaM-Full_torc display_gen

AMaLGaM-Full-Free: AMaLGaM-Full-Free.c
	$(CC) $(CFLAGS) -o AMaLGaM-Full-Free AMaLGaM-Full-Free.c $(LIBS)

AMaLGaM-Full-Free_torc: AMaLGaM-Full-Free_torc.c
	$(TCC) $(TCFLAGS) -o AMaLGaM-Full-Free_torc AMaLGaM-Full-Free_torc.c $(TLIBS)

AMaLGaM-Full_torc: AMaLGaM-Full_torc.c
	$(TCC) $(TCFLAGS) -o AMaLGaM-Full_torc AMaLGaM-Full_torc.c $(TLIBS)

display_gen: display_gen.c gnuplot_i.c
	$(CC) $(CFLAGS) -o display_gen display_gen.c gnuplot_i.c $(LDFLAGS)

clean:
	rm -f AMaLGaM-Full-Free AMaLGaM-Full-Free_torc AMaLGaM-Full_torc display_gen

clear:
	rm -f *.dat

cleanall:
	rm -f *.dat

# ./AMaLGaM-Full-Free_torc -w -v 12 2 -6 6 0                    10000 0    0.000001 1
# ./AMaLGaM-Full_torc            12 2 -6 6 0 0.35 25 10 0.9 1.0 10000 0 10 0.000001
