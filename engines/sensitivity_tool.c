/*
 *
 *  sensitivity_tool.c
 *  Pi4U
 *
 *  Created by Line Kulakova on 23/8/16.
 *  Copyright 2016 ETH Zurich. All rights reserved.
 *
 */

#include "engine_tmcmc.h"
#include "fitfun.c"

#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <torc.h>
#include <math.h>

void taskfun(double *x, int *pn, double *res, int *info)
{
    int n = *pn;
	*res = fitfun(x, n, (void *)NULL, info);
    return;
}


#define PROBDIM	2
#define NRUNS 5 // how many steps to do each side from the reference point

int main(int argc, char *argv[])
{
	int i, j, k;

	torc_register_task(taskfun);

	torc_init(argc, argv, MODE_MS);


	char ifilename[80];
	strcpy(ifilename, "reference_point.txt");

	if (argc > 1) strcpy(ifilename, argv[1]);

	FILE *fp = fopen(ifilename, "r");
	if (fp == NULL) {
		printf("Input file %s does not exist!\n", ifilename);
		torc_finalize();
		return 1;
	}

	double ref_t[PROBDIM], ref, res_t[PROBDIM][2*NRUNS], res[PROBDIM][2*NRUNS];
	for (i=0; i<PROBDIM; ++i)
		fscanf(fp, "%lf", &ref_t[i]);
	fclose(fp);

	printf("REFERENCE: ");
	for (i=0; i<PROBDIM; ++i)
		printf("%.2lf ", ref_t[i]);
	printf("\n");

	int n = PROBDIM;

	// reference case
	int info[4];
	info[0] = 0; info[1] = 0; info[2] = 0; info[3] = -1;
	torc_create(-1, taskfun, 4,
		PROBDIM, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_INT, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_RES,
		4, MPI_INT, CALL_BY_COP,
		ref_t, &n, &ref, info);

	// run simulations on a mesh for each parameter
	for (i=0; i<PROBDIM; ++i)
	{
		for (j=1; j<=NRUNS; ++j)
		{
			res_t[i][NRUNS-j] = (1 - 0.01*j)*ref_t[i]; // subtract j percent

			double t[PROBDIM];
			for (k=0; k<PROBDIM; ++k)
				t[k] = ref_t[k];
			t[i] = res_t[i][NRUNS-j];

			printf("t: ");
			int p;
			for (p=0; p<PROBDIM; ++p)
				printf("%.4lf ", t[p]);
			printf("\n");

			int info[4];
			info[0] = 0; info[1] = 0; info[2] = 0;
			info[3] = i*2*NRUNS + (NRUNS-j);

			torc_create(-1, taskfun, 4,
				PROBDIM, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_RES,
				4, MPI_INT, CALL_BY_COP,
				t, &n, &res[i][NRUNS-j], info);
		}

		for (j=1; j<=NRUNS; ++j)
		{
			res_t[i][NRUNS+j-1] = (1 + 0.01*j)*ref_t[i]; // add j percent

			double t[PROBDIM];
			for (k=0; k<PROBDIM; ++k)
				t[k] = ref_t[k];
			t[i] = res_t[i][NRUNS+j-1];

			printf("t: ");
			int p;
			for (p=0; p<PROBDIM; ++p)
				printf("%.4lf ", t[p]);
			printf("\n");

			int info[4];
			info[0] = 0; info[1] = 0; info[2] = 0;
			info[3] = i*2*NRUNS + (NRUNS+j-1);

			torc_create(-1, taskfun, 4,
				PROBDIM, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_RES,
				4, MPI_INT, CALL_BY_COP,
				t, &n, &res[i][NRUNS+j-1], info);
		}
	}
	torc_waitall();

	for (i=0; i<PROBDIM; ++i) // loop over parameters
	{
		for (j=0; j<2*NRUNS; ++j)
		{
			printf("RESULT ( ");
			for (k=0; k<i; ++k) printf("%.16lf ", ref_t[k]);
			printf("%.16lf ", res_t[i][j]);
			for (k=i+1; k<PROBDIM; ++k) printf("%.16lf ", ref_t[k]);
			printf(") = ");
			printf("%.16lf\n", res[i][j]);
		}
	}

	torc_finalize();

	return 0;
}

