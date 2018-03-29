#ifndef _LOGLIKE_THETA_H_
#define _LOGLIKE_THETA_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void loglike_theta_initialize();
void loglike_theta_finalize();
double loglike_theta(double *x, int n, void *output, int *winfo);

double loglike_(double *x, int n, void *output, int *winfo);

#endif
