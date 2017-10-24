/*
 *  display_gen.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "gnuplot_i.h"

int display = 1;

gnuplot_ctrl * g = NULL;

char *FILENAME;
// "curgen"
// "samples"
// "seeds"

void display_prop(int nsamples)
{
        FILE *fp;
        char fname[256];

	sprintf(fname, "%s", FILENAME);
        fp = fopen(fname, "r");
        if (fp == NULL) {
                printf("No file %s\n", fname);
		return;
        }
	fclose(fp);

	if (!display) {
		gnuplot_cmd(g, "set terminal png");
		gnuplot_cmd(g, "set output \"%s.png\"", fname);
	}


	{
	char title_str[64];
	sprintf(title_str, "\"%d_%d_%d\"", 0,1,1);
	gnuplot_cmd(g, "set term x11 %d persist title %s", 0, title_str);

	gnuplot_cmd(g, "set nokey");
	int i;
	for (i = 1; i <= nsamples; i++) {
		char cmd[256];
		if (i == 1)
			sprintf(cmd, "plot \"%s\" using %d pt 2", FILENAME, i);
		else
			sprintf(cmd, "replot \"%s\" using %d pt 2", FILENAME, i);
		gnuplot_cmd(g, cmd); //"splot \"%s\" %s with points pt 7 palette", fname, using_str);
	}
	}

	sleep(1);
}


int main(int argc, char *argv[])
{
	int nsamples = 0;

	if (argc != 3) {
		printf("usage: %s <filename> <#samples>\n", argv[0]);
		exit(1);
	}

	FILENAME = argv[1];
	nsamples = atoi(argv[2]);

	g = gnuplot_init();
	display_prop(nsamples);
	gnuplot_close(g);

	return 0;
}
