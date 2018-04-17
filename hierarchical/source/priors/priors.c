#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "priors.h"
#include "myrand.h"

#define BLANKS " \t"

//=============================================================================
//
//
double eval_density(Density d, double x){
	return d.f(x,d.par);
}



double eval_log_density(Density d, double x){
	return d.lf(x,d.par);
}


double eval_random( Density d ){
	return d.r( d.par );
}


double print_density( Density d ){
	
	if( strcmp(d.name,"uniform")==0 )
		printf("%s: %lf  -  %lf \n", d.name, d.par[0],d.par[1]);
	else if( strcmp(d.name,"normal")==0 )
		printf("%s: %lf  -  %lf \n", d.name, d.par[0],d.par[1]);
	else if( strcmp(d.name,"exponential")==0 )
		printf("%s: %lf \n", d.name, d.par[0]);
	else if( strcmp(d.name,"gamma")==0 )
		printf("%s: %lf  -  %lf \n", d.name, d.par[0],d.par[1]);

}


double print_priors( Density *d, int N){
	printf("==============================\n");
	printf("===  Prior Distribution    ===\n");
	printf("==============================\n");

	for(int i=0; i<N; i++)
		print_density( d[i] );

	printf("==============================\n");
}




double prior_pdf( Density *d, int N, double x){

	double res=1.;
	
	for(int i=0; i<N; i++)
		res *= eval_density( d[i], x);
}



double prior_log_pdf( Density *d, int N, double x){

	double res=0.;
	
	for(int i=0; i<N; i++)
		res += eval_log_density( d[i], x);
}




//=============================================================================
//
//
void read_priors( Density **p_priors, int *p_N ){

	FILE *fp = fopen("priors.par","r");
	if(fp==NULL){
		printf("\npriors.par could not be opened. Exit...\n");
		exit(EXIT_FAILURE);
	}

	
	int N=-1;	


	size_t bufsize = 1024, nchar=1;
	char * buffer = (char *)malloc(bufsize * sizeof(char));	
	int cnt = 0;

	Density *prior=NULL;
	
	while(  -1 != getline( &buffer, &bufsize, fp) ){
			
		char * pch = strtok (buffer, BLANKS );
  		while( pch != NULL ){
			
			if( pch[0]=='#' || pch[0]=='\n' ) //go to the next line. 
				break;
			

			if( strcmp(pch,"N")==0 ){
    			pch = strtok (NULL, BLANKS );
				N = atoi(pch);
				check_n(N);
				prior = (Density *)malloc(N*sizeof(Density));
				break;

			}


			if( strcmp(pch,"uniform")==0 || strcmp(pch,"uni")==0 ){
				check_n(N);
				prior[cnt].f  = &uniform_pdf;
				prior[cnt].lf = &uniform_log_pdf;
				prior[cnt].r  = &uniform_rnd;


				prior[cnt].par = (double*)malloc(2*sizeof(double));
				pch = strtok (NULL, BLANKS );
				prior[cnt].par[0] = atof(pch);
				pch = strtok (NULL, BLANKS );
				prior[cnt].par[1] = atof(pch);

				strcpy( prior[cnt].name, "uniform" );
				cnt++;
				break;
			}


			if( strcmp(pch,"normal")==0 || strcmp(pch,"gaussian")==0 ){
				check_n(N);
				prior[cnt].f  = &normal_pdf;
				prior[cnt].lf = &normal_log_pdf;
				prior[cnt].r  = &normal_rnd;


				prior[cnt].par = (double*)malloc(2*sizeof(double));
				pch = strtok (NULL, BLANKS );
				prior[cnt].par[0] = atof(pch);
				pch = strtok (NULL, BLANKS );
				prior[cnt].par[1] = atof(pch);

				strcpy( prior[cnt].name, "normal" );
				cnt++;
				break;
			}
			
			if( strcmp(pch, "exp")==0 || strcmp(pch,"exponential")==0  ){
				check_n(N);
				prior[cnt].f  = &exp_pdf;
				prior[cnt].lf = &exp_log_pdf;
				prior[cnt].r  = &exp_rnd;


				prior[cnt].par = (double*)malloc(1*sizeof(double));
				pch = strtok (NULL, BLANKS );
				prior[cnt].par[0] = atof(pch);
				
				strcpy( prior[cnt].name, "exponential" );
				cnt++;
				break;
			}

			if( strcmp(pch,"gamma")==0 ){
				check_n(N);
				prior[cnt].f  = &gamma_pdf;
				prior[cnt].lf = &gamma_log_pdf;
				prior[cnt].r  = &gamma_rnd;


				prior[cnt].par = (double*)malloc(2*sizeof(double));
				pch = strtok (NULL, BLANKS );
				prior[cnt].par[0] = atof(pch);
				pch = strtok (NULL, BLANKS );
				prior[cnt].par[1] = atof(pch);

				strcpy( prior[cnt].name, "gamma" );
				cnt++;
				break;
			}

			puts(buffer);
			printf("\nSomething went wrong while reading the priors parameter file. Exit...\n");
			exit(EXIT_FAILURE);

		}
	}

	fclose(fp);

	*p_priors = prior;
	*p_N = N;
}





//=============================================================================
//
//
void check_n( int N ){
	if(N<=0){
		puts("Check the prios.par file. Error with reading N. Exit...\n");
		exit(EXIT_FAILURE);
	}
}


