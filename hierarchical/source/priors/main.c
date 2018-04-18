#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "priors.h"
#include "myrand.h"


void test_suite( Density *priors, int N);
double mean_std( double *data, int N, double *mean, double *std );




int main ()
{

	int N;
	Density *priors;

	read_priors( &priors, &N );

	gsl_rand_init( 123456 );

	print_priors( priors, N);

	test_suite( priors, N);
  
	return 0;
}











void test_suite( Density *p, int N){



	for(int i=0; i<N; i++){
		
		printf("%s: \n", p[i].name );
		
		double v1 = eval_density(p[i], 1.0);
		double v2 = eval_log_density(p[i], 1.0);
		
		printf("	%lf - %lf \n", v1, exp(v2));
	
	}

	printf("-----------\n");
	printf("Prior: \n");
	
	double x[4]={1,1,1,1};
	double v1 = prior_pdf( p, N, x);
	double v2 = prior_log_pdf( p, N, x);
	printf("	%lf - %lf \n", v1, exp(v2));

	printf("-----------\n");

	int Ns = 1e6;
	double rnd[Ns], mean, std, meanx, stdx;
	
	for(int i=0; i<N; i++){

		char s[64];
		sprintf(s, "rnd_%02d.txt",i);
		FILE *fp = fopen(s,"w");
		
		for(int j=0; j<Ns; j++){
			rnd[j] = eval_random( p[i]);
			fprintf(fp,"%.12lf \n", rnd[j] );
		}

		mean_std(rnd, Ns, &mean, &std);

		if(strcmp(p[i].name,"normal")==0){
			meanx = p[i].par[0];
			stdx  = p[i].par[1];
		}
		else if(strcmp(p[i].name,"uniform")==0){
			meanx = (p[i].par[0]+p[i].par[1]) / 2;
			stdx  = (p[i].par[1]-p[i].par[0]) / sqrt(12);
		}
		else if(strcmp(p[i].name,"exponential")==0){
			meanx = p[i].par[0];
			stdx  = p[i].par[0];
		}
		else if(strcmp(p[i].name,"gamma")==0){
			meanx = p[i].par[0] * p[i].par[1];
			stdx  = sqrt(p[i].par[0]) * p[i].par[1];
		}


		printf( "%s: \n", p[i].name, mean, std );
		printf( "	exact:    %lf  ---  %lf \n", meanx, stdx );
		printf( "	computed: %lf  ---  %lf \n", mean, std );

	


		fclose(fp);


	}	


}









double mean_std( double *data, int N, double *mean, double *std ){
    
	*mean = 0.;
	*std  = 0.;

    for( int i=0; i<N; i++)
        *mean += data[i];

    *mean = (*mean)/N;

    for( int i=0; i<N; i++)
        *std += pow( data[i] - (*mean) , 2);

	*std = sqrt( *std / (N-1) );
}












