/*
 *  engine_tmcmc.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <ftw.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include "fitfun.c"
#include <mpi.h>
#include <torc.h>


void taskfun(double *TP, int *pn, int info[4], double *res)
{
	int i, n = *pn;
	*res = fitfun(TP, n, 0, info);
}

#define PROBDIM 2

int main(int argc, char *argv[])
{
	int i;
	double TP[PROBDIM];

	torc_register_task(taskfun);
	torc_init(argc, argv, MODE_MS);

	srand48(10);

	int n = PROBDIM;
	int t;

	double base[PROBDIM] = {0.52576, 0.0331156};	// one initial point

	#define NTASKS	(6)
	double res[NTASKS];
	int info[4];
	for (t = 0; t < NTASKS; t++) {
		for (i = 0; i < PROBDIM; i++) {
			double d = drand48();
			if (d > 0.5)
				TP[i] = base[i]*(1 + d*0.01);
			else
				TP[i] = base[i]*(1 - d*0.01);			

		}
		info[0] = 0; info[1] = 0; info[2] = 0; info[3] = t;

		torc_create(-1, taskfun, 4,
			PROBDIM, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_INT, CALL_BY_COP,
			4, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			TP, &n, info, &res[t]);
	}
	torc_waitall();

	for (t = 0; t < NTASKS; t++)
		printf("res[%d] = %lf\n", t, res[t]);

	torc_finalize();

	return 0;
}

