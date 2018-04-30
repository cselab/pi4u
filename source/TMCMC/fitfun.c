#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


static gsl_vector *mu;
static gsl_matrix *sigma;


void fitfun_initialize(char *s){
	
	int N = atoi(s);

	mu = gsl_vector_calloc(N);
	for( int i=0; i<N; i++)
		gsl_vector_set( mu, i, (double)i );


	//gsl_vector_fprintf( stdout, mu, "%f");


	sigma = gsl_matrix_calloc( N,N );

	double D0;
	for( int i=0; i<N; i++){
		D0 = (i%2)*0.5 + 0.02;
		gsl_matrix_set( sigma, i, i, D0 );
	}

	for( int i=1; i<N; i++){
		D0 = (2*(i%2) - 1) * 0.05; //change the sign
		gsl_matrix_set( sigma, i, i-1, D0 );
	}


	//printf("----->\n");
	//gsl_matrix_fprintf( stdout, sigma, "%f");

	gsl_linalg_cholesky_decomp( sigma );


}



void fitfun_finalize(){

	gsl_vector_free(mu);
	gsl_matrix_free(sigma);

}









double fitfun(double *x, int N, void *output, int *info){


	gsl_vector *xg 	= gsl_vector_calloc(N);
	gsl_vector *work = gsl_vector_calloc(N);

	double result;

	for( int i=0; i<N; i++)
		gsl_vector_set( xg, i, x[i] );

	gsl_ran_multivariate_gaussian_log_pdf( xg, mu, sigma, &result, work);

	gsl_vector_free(xg);
	gsl_vector_free(work);

	return result;

}




/*
double fitfun(double *x, int N, void *output, int *info){


	double  v1=0,v2=0;

	for(int i=0; i<N; i++){
		v1 += log(gsl_ran_gaussian_pdf( x[i],1.0 ));
		v2 += log(gsl_ran_gaussian_pdf( x[i]-6,2.0 ));
	}

	return log( 0.5*exp(v1) + 0.5*exp(v2)  )  ;

}
*/

