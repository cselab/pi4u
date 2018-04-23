#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_randist.h>



void fitfun_initialize(char *s){

}



void fitfun_finalize(){

}



double fitfun(double *x, int N, void *output, int *info){


	double  v1=0,v2=0;

	for(int i=0; i<N; i++){
		v1 += log(gsl_ran_gaussian_pdf( x[i],1.0 ));
		v2 += log(gsl_ran_gaussian_pdf( x[i]-6,2.0 ));
	}

	return log( 0.5*exp(v1) + 0.5*exp(v2)  )  ;

}

