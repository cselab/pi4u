#include <stdio.h>
#include <string.h>
#include <Rinternals.h>
#include <Rembedded.h>

char fitfun_name[256];

void fitfun_initialize(char *fname)
{
	strcpy(fitfun_name, fname);
}

void fitfun_finalize()
{
}


void source(const char *name)
{
	SEXP e;

	PROTECT(e = lang2(install("source"), mkString(name)));
	R_tryEval(e, R_GlobalEnv, NULL);
	UNPROTECT(1);
}

/**
 * Wrapper for R function fitfun, defined in fitfun.R.
 */
double R_fitfun(double a[], int alen)
{
	// Allocate an R vector and copy the C array into it.
	SEXP arg;
	PROTECT(arg = allocVector(REALSXP, alen));
	memcpy(REAL(arg), a, alen * sizeof(double));

	// Setup a call to the R function
	SEXP fitfun_call;
	//PROTECT(fitfun_call = lang2(install("fitfun"), arg));
	PROTECT(fitfun_call = lang2(install(fitfun_name), arg));

	// Execute the function
	int errorOccurred;
	SEXP ret = R_tryEval(fitfun_call, R_GlobalEnv, &errorOccurred);

	double res;

	if (!errorOccurred)
	{
//		printf("R returned: ");
		double *val = REAL(ret);
//		for (int i = 0; i < LENGTH(ret); i++)
//			printf("%0.1f, ", val[i]);
//		printf("\n");

		res = val[0];
	}
	else
	{
		printf("Error occurred calling R\n");
		res = -1e12;
	}

	UNPROTECT(2);

	return res;
}


double fitfun(double *x, long N, void *output, int *info)
{
	char fitfun_file[256];

	sprintf(fitfun_file,"%s.R", fitfun_name);
//	source("fitfun.R");
	source(fitfun_file);
	double r = R_fitfun(x, N);
	return r;
}

#include "engine_tmcmc.h"

void tmcmcR(double *res, int *tmcmc_info, int *pNth, int *pMaxStages, int *pPopSize, double *lb, double *ub)
{
	int Nth = *pNth;
	int MaxStages = *pMaxStages;
	int PopSize = *pPopSize;

	int i;

	for (i = 0; i < 4; i++) printf("tmcmc_info[%d]=%d\n", i, tmcmc_info[i]);
	printf("Nth = %d\n", Nth);
	printf("MaxStages = %d\n", MaxStages);
	printf("PopSize = %d\n", PopSize);
	for (i = 0; i < Nth; i++)
	{
		printf("lb[%d]=%f\n", i, lb[i]);
	}
	for (i = 0; i < Nth; i++)
	{
		printf("ub[%d]=%f\n", i, ub[i]);
	}

  	tmcmc(res, tmcmc_info, Nth, MaxStages, PopSize, lb, ub);
	return;
}

void tmcmcR_initialize(char **fname)
{
	int i = tmcmc_initialize(*fname);
}

void tmcmcR_finalize()
{
	tmcmc_finalize();
}
