all: normaltest bananatest generaltest

normaltest: normaltest.c dramrun.c
	gcc -Wall -std=c99 -o normaltest normaltest.c -lgsl -lm
#	gcc -Wall -DDEBUG=1 -std=c99 -o normaltest normaltest.c -lgsl -lm


bananatest: bananatest.c dramrun.c
	gcc -Wall -std=c99 -o bananatest bananatest.c -lgsl -lm

generaltest: generaltest.c dramrun.c fitfun.c
	gcc -Wall -std=c99 -o generaltest generaltest.c -lgsl -lm

clean:
	rm -f normaltest bananatest generaltest

clear:
	rm -f *test_chain.txt
