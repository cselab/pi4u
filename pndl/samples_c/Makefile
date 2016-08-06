all: psample0.exe psample1.exe

MPIF90 = mpif90 -Wall
MPICC = mpicc -Wall

psample0.exe: psample0.o ../src/libpndl.a 
	$(MPIF90) -o psample0.exe psample0.o  -L../src  -lpndl `torc_libs` 

psample0.o: psample0.c
	$(MPICC) `torc_cflags` -I../src  -O3 -c  psample0.c -o psample0.o

psample1.exe: psample1.o ../src/libpndl.a 
	$(MPIF90) -o psample1.exe psample1.o  -L../src  -lpndl `torc_libs` 

psample1.o: psample1.c
	$(MPICC) `torc_cflags` -I../src -O3 -c  psample1.c -o psample1.o

clean:
	rm -f *.o *.exe
