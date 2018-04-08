#include "../TMCMC/fitfun.h"
#include "loglike_theta_fast.h"

void fitfun_initialize() {

	loglike_theta_fast_initialize();

}



void fitfun_finalize() {

	loglike_theta_fast_finalize();

}




double fitfun(double *x, int n, void *output, int *info) {

	return loglike_theta_fast( x, n, output );

}
