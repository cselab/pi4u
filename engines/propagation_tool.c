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
#include <stdlib.h>
#include <mpi.h>
#include <torc.h>
#include <math.h>

void taskfun(double *x, int *pn, double *res, int *info)
{
        int n = *pn;
	int gen, chain, step, task;
	gen = info[0]; chain = info[1]; step = info[2]; task = info[3];
	printf("executing task (%d,%d,%d,%d,*)\n", gen, chain, step, task);

#if 1
	*res = fitfun(x, n, (void *)NULL, info);

#else
	/* this is a very simple example */
	static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
	static FILE *fp = NULL;

	pthread_mutex_lock(&m);
	if (fp == NULL) fp = fopen("propagation.txt", "w");
	pthread_mutex_unlock(&m);

        double obs[12]; /* e.g. assuming 12 observations */
        *res = fitfun(x, n, (void *)obs, info);

	pthread_mutex_lock(&m);
        int i;
	for (i = 0; i < 12; i++) fprintf(fp, "%f ", obs[i]);
        fprintf(fp, "\n");
	pthread_mutex_unlock(&m);

#endif

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
	char line[2048];
	double TP[MAXPOINTS][PROBDIM], ref[MAXPOINTS], res[MAXPOINTS];
	int t = 0;
	while (fgets(line, 2048, fp)!= NULL)
	{
		printf("line %d: %s", t, line);
		char *p;
		p = strtok(line, " \t\n");
		for (i = 0; i < PROBDIM; i++) {
			if (!p) break;
			TP[t][i] = atof(p);
			p = strtok(NULL, " \t\n");
		}
		if (p) p = strtok(NULL, " \t\n");
		if (p) ref[t] = atof(p);

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

