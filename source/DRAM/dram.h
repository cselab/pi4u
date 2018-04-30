/*
 *  dram.h
 *  Pi4U
 *
 *  Copyright 2018 ETH Zurich. All rights reserved.
 *
 */

#ifndef _DRAM_H_
#define _DRAM_H_

#include "../priors/priors.h"


/* Variables related to the inferred parameters */
struct _params
{
	double 	*par0;		/* Initial parameter vector */
	double 	sigma2;		/* Initial/prior value for the Gaussian error variance */
	double 	n0;			/* Precision of sigma2 as imaginative observations. if n0<0, no sigma2 update */
	int 	n;			/* Nnumber of actual observations (for sigma2 update) */

	double 	*lbounds;	/* Lower bounds of parameter values */
	double 	*ubounds;	/* Upper bounds of parameter values */

	Density *prior;

	char 	filename[256];	/* Output file for the chain */

} params;



/* Variables related to the options of DRAM */
struct _options
{
	int 	Npar;    	/* problem dimensionality (number of parameters) */
    int 	Nsim;    	/* Length of the chain */

    double 	DRscale;    /* DR shrink factor, if zero, no DR */

	int 	AMinterv;	/* how often to adapt, if zero, no adaptation */	// 20 //= 500;
	double 	AMscale;	/* Scale for adapting the proposal = 2.4/sqrt(Npar) */
	double 	AMepsilon;	/* Small constant used in the AM covariance update */

	double 	*qcov;     	/* Proposal covariance matrix */

	int 	printfreq;	/* Time interval for printing info on screen */
	int 	verbose;

} options;



void dram_init();
void dram();
void dram_finalize();



#endif
