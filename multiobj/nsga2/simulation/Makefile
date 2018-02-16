CC = gcc

all: simcode install

simcode: simcode.c
	$(CC) -o simcode simcode.c -lm

install:
	mkdir -p ../bin
	cp -f simcode ../bin
	cp -f run.sh ../bin

clean:
	rm -f simcode params.txt fitness.txt constr_viol.txt

