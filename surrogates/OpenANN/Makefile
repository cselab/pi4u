all:	mytest

dim?=2

CFLAGS=-Wall -O3
#CFLAGS+=-g -ggdb
#CFLAGS+=-DVERBOSE=1

nnet.o: nnet.cpp
	g++ $(CFLAGS) nnet.cpp -c -o nnet.o `pkg-config --cflags openann`

mytest: mytest.c nnet.o
	gcc $(CFLAGS) -DDIM=$(dim) mytest.c -o mytest nnet.o `pkg-config --cflags --libs openann` -lalglib -lstdc++

clean:
	rm -f mytest nnet.o myNN.txt
