#include <fitfun.h>
#include "loglike_theta.h"


void fitfun_initialize(char *s) {

}



void fitfun_finalize() {

}




double fitfun(double *x, int n, void *output, int *info) {

	return loglike_theta( x, n, output, info );

}
