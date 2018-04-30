#ifndef _LOGLIKE_PSI_H_
#define _LOGLIKE_PSI_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



void 	loglike_psi_initialize();
void 	loglike_psi_finalize();
double 	loglike_psi(double *x, int n, void *output, int *info);

double log_priorHB(double *x, double *psi, int n);


#endif
