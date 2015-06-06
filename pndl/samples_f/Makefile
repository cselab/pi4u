all: psample3s.exe psample3a.exe 

MPIF90 = mpif90
MPICC = mpicc

psample3s.exe: psample3s.o ../src/libpndl.a 
	$(MPIF90) -o psample3s.exe psample3s.o  -L../src  -lpndl `torc_libs`

psample3s.o: psample3s.f
	$(MPIF90) `torc_cflags`  -O3 -c  psample3s.f -o psample3s.o

psample3a.exe: psample3a.o ../src/libpndl.a 
	$(MPIF90) -o psample3a.exe psample3a.o  -L../src  -lpndl `torc_libs`

psample3a.o: psample3a.f
	$(MPIF90) `torc_cflags`  -O3 -c  psample3a.f -o psample3a.o

clean:
	rm -f *.o *.exe
