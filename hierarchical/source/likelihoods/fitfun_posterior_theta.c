#include <fitfun.h>
#include "loglike_posterior_theta.h"



void fitfun_initialize(char *s) {

	loglike_posterior_theta_initialize();

}



void fitfun_finalize() {

	loglike_posterior_theta_finalize();

}




double fitfun(double *x, int n, void *output, int *info) {

	return loglike_posterior_theta( x, n, output, info );

}
