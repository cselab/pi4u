CC=mpicc
CFLAGS=`torc_cflags` -O3 -D_XOPEN_SOURCE=700 -D_BSD_SOURCE

LIBS = `torc_libs`

all: driver_rosen.exe simcode

driver_rosen.exe: driver_rosen.o 
	$(CC) $(CFLAGS) -o driver_rosen.exe driver_rosen.o $(LIBS)

driver_rosen.o: driver_rosen.c 
	mpicc $(CFLAGS)  -c  driver_rosen.c 

simcode:
	cd bin; make


clean:
	cd bin; make clean
	rm -f *.o *.exe
