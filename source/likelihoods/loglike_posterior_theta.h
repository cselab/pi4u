#ifndef _LOGLIKE_POSTERIOR_THETA_FAST_H_
#define _LOGLIKE_POSTERIOR_THETA_FAST_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void 	loglike_posterior_theta_initialize();
void 	loglike_posterior_theta_finalize();
double 	loglike_posterior_theta(double *x, int n, void *output, int *info );

double 	loglike_(double *x, int n, void *output, int *info );

#endif
