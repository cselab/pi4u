/*
 *  propagation_tool.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 29/6/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#include "fitfun.c" 
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <torc.h>
#include <math.h>

void taskfun(double *x, int *pn, double *res, int *info)
{
        int n = *pn;
	int gen, chain, step, task;
	gen = info[0]; chain = info[1]; step = info[2]; task = info[3];
	printf("executing task (%d,%d,%d,%d,*)\n", gen, chain, step, task);

	*res = fitfun(x, n, (void *)NULL, info);
        return;
}


#define PROBDIM	2

int main(int argc, char *argv[])
{
	int i;

	torc_register_task(taskfun);

	torc_init(argc, argv, MODE_MS);

	int n = PROBDIM;

	char ifilename[80];
	strcpy(ifilename, "points.txt");

	if (argc > 1) strcpy(ifilename, argv[1]);

	FILE *fp = fopen(ifilename, "r");
	if (fp == NULL) {
		printf("Input file %s does not exist!\n", ifilename);
		torc_finalize();
		return 1;
	}

/*
	int nlines=0;
	char line[256];	
	while (fgets(line, 256, fp)!= NULL) {
		nlines++;
	}
*/

	#define MAXPOINTS 64

	char line[256];
	double TP[MAXPOINTS][PROBDIM], ref[MAXPOINTS], res[MAXPOINTS];
	int t = 0;
	while (fgets(line, 256, fp)!= NULL)
	{
		sscanf(line, "%lf %lf %lf", &TP[t][0], &TP[t][1], &ref[t]);
		printf("line %d: %s", t, line);
		t++;
		if (t == MAXPOINTS) break;
	}
	fclose(fp);

	for (i = 0; i < t; i++) {
		int info[4];
		info[0] = 0; info[1] = 0; info[2] = 0; info[3] = i;

		torc_create(-1, taskfun, 4,
			PROBDIM, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			4, MPI_INT, CALL_BY_COP,
			TP[i], &n, &res[i], info);
	}
	torc_waitall();

	for (i = 0; i < t; i++) {
		printf("RESULT %03d: %10.4f %10.4f %10.4e %10.4e %10.4lf\n", i, TP[i][0], TP[i][1], ref[i], res[i], fabs(ref[i]-res[i]));
	}

	torc_finalize();

	return 0;
}

