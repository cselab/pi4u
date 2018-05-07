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



void delete_density( Density *d ){
	free(d->par);
}


void delete_prior( Density *d, int N){
    
    if( d==NULL )
        return;

	for( int i=0; i<N; i++ )
		delete_density( &d[i] );
	
	free(d);
}



double eval_density(Density d, double x){
	return d.f(x,d.par);
}



double eval_log_density(Density d, double x){
	return d.lf(x,d.par);
}


double eval_random( Density d ){
	return d.r( d.par );
}


void print_density( Density d ){
	
	if( strcmp(d.name,"uniform")==0 )
		printf("%s: %lf  -  %lf \n", d.name, d.par[0],d.par[1]);
	else if( strcmp(d.name,"normal")==0 )
		printf("%s: %lf  -  %lf \n", d.name, d.par[0],d.par[1]);
	else if( strcmp(d.name,"exponential")==0 )
		printf("%s: %lf \n", d.name, d.par[0]);
	else if( strcmp(d.name,"gamma")==0 )
		printf("%s: %lf  -  %lf \n", d.name, d.par[0],d.par[1]);

}


void print_priors( Density *d, int N){
	printf("==============================\n");
	printf("===  Prior Distribution    ===\n");
	printf("==============================\n");

	for(int i=0; i<N; i++)
		print_density( d[i] );

	printf("==============================\n");
}




double prior_pdf( Density *d, int N, double *x){

	double res=1.;
	
	for(int i=0; i<N; i++)
		res *= eval_density( d[i], x[i] );
	
	return res;
}



double prior_log_pdf( Density *d, int N, double *x){

	double res=0.;
	
	for(int i=0; i<N; i++)
		res += eval_log_density( d[i], x[i] );
	
	return res;
}




//=============================================================================
//
//
void read_priors( const char *file, Density **p_priors, int *p_N ){

	FILE *fp = fopen( file,"r");
	if(fp==NULL){
		printf("\n%s could not be opened. Exit...\n", file );
		exit(EXIT_FAILURE);
	}

	
	int N=-1;	


	size_t bufsize = 1024;
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
				prior[cnt].npar = 2;
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
				prior[cnt].npar = 2;
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
				prior[cnt].npar = 1;
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
				prior[cnt].npar = 2;
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
	free(buffer);

	*p_priors = prior;
	*p_N = N;
}










//=============================================================================
//
//
void new_prior_from_prior( Density **new_prior, Density *from_prior, int Npr ){


	Density *prior = (Density *) malloc( Npr*sizeof(Density) );

	for( int i=0; i<Npr; i++ ){

		prior[i].f  = from_prior[i].f;
		prior[i].lf = from_prior[i].lf;
		prior[i].r  = from_prior[i].r;

		prior[i].npar = from_prior[i].npar;
		prior[i].par  = (double*) malloc( sizeof(double)*prior[i].npar );

		for( int j=0; j<prior[i].npar; j++)
			prior[i].par[j] = from_prior[i].par[j];

		strcpy( prior[i].name, from_prior[i].name );

	}

	*new_prior = prior;
}












//=============================================================================
//
//
//	XXX: danger: the loop may go outside of psi. no length check. chould we pass the length of
//	psi as argument? 
void reassign_prior( Density *p, int Np, double *psi ){

	int cnt = 0;
	for( int i=0; i<Np; i++)
		for( int j=0; j<p[i].npar; j++)
				p[i].par[j] = psi[cnt++];

}






//=============================================================================
//
//
void check_n( int N ){
	if(N<=0){
		puts("Check the priors parameter file. Error with reading N. Exit...\n");
		exit(EXIT_FAILURE);
	}
}


