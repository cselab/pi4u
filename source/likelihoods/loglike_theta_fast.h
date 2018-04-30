#ifndef _LOGLIKE_THETA_FAST_H_
#define _LOGLIKE_THETA_FAST_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void 	loglike_theta_fast_initialize();
void 	loglike_theta_fast_finalize();
double 	loglike_theta_fast(double *x, int n, void *output );

double 	loglike_(double *x, int n, void *output, int *info );


#endif
