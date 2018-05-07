#define _GNU_SOURCE

#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#include <math.h>

#include <fitfun.h>
#include "loglike_posterior_theta.h"
#include "../priors/priors.h"


#define BLANKS " \t\n="
#define BUFLEN 1024


typedef struct database_s{
	char theta_folder[128];
	char theta_file[128];
	char logev_file[128];
	
	char psi_folder[128];
	char psi_file[128];

	int  Ntheta;
	int  Npsi;

	int Dtheta;
	int Dpsi;

}database;


static database db;

static double logev_theta;
static double *denom;



static Density **priors;
static int Npr;


//=============================================================================
//
//
// count the number of lines in file fp  and check if less than N
static int enough_lines(FILE *fp, int N) {
    char ch;
    int lines = 0;

	rewind(fp);

    while (!feof(fp)){
            ch = fgetc(fp);
            if (ch == '\n')
                lines++;
    }
    rewind(fp);

    return lines >= N;
}




//=============================================================================
//
//
void read_db_file( const char *file  ){

	FILE *fp = fopen( file, "r" );
	if(fp==NULL){
		printf("\n%s could not be opened. Exit...\n", file );
		exit(EXIT_FAILURE);
	}

	
	size_t bufsize = 1024;
	char * buffer = (char *)malloc(bufsize * sizeof(char));	

	while(  -1 != getline( &buffer, &bufsize, fp) ){
			
		char * pch = strtok (buffer, BLANKS );
  		while( pch != NULL ){
			
			if( pch[0]=='#' || pch[0]=='\n' ) //go to the next line. 
				break;
			
			if( strcmp(pch,"theta_folder")==0 ){
    			pch = strtok (NULL, BLANKS );
				strcpy( db.theta_folder, pch );
				break;
			}

			if( strcmp(pch,"psi_folder")==0 ){
    			pch = strtok (NULL, BLANKS );
				strcpy( db.psi_folder, pch );
				break;
			}

			if( strcmp(pch,"theta_file")==0 ){
    			pch = strtok (NULL, BLANKS );
				strcpy( db.theta_file, pch );
				break;
			}


			if( strcmp(pch,"logev_file")==0 ){
    			pch = strtok (NULL, BLANKS );
				strcpy( db.logev_file, pch );
				break;
			}

			if( strcmp(pch,"psi_file")==0 ){
    			pch = strtok (NULL, BLANKS );
				strcpy( db.psi_file, pch );
				break;
			}

			if( strcmp(pch,"Ntheta")==0 ){
    			pch = strtok (NULL, BLANKS );
				db.Ntheta = atoi( pch );
				break;
			}

			if( strcmp(pch,"Npsi")==0 ){
    			pch = strtok (NULL, BLANKS );
				db.Npsi = atoi( pch );
				break;
			}

			if( strcmp(pch,"Dtheta")==0 ){
    			pch = strtok (NULL, BLANKS );
				db.Dtheta = atoi( pch );
				break;
			}

			if( strcmp(pch,"Dpsi")==0 ){
    			pch = strtok (NULL, BLANKS );
				db.Dpsi = atoi( pch );
				break;
			}


			puts(buffer);
			printf("\nSomething went wrong while reading the theta database file %s. Exit...\n",file);
			exit(EXIT_FAILURE);

		}
	}

	fclose(fp);
	
	free(buffer);
}





//=============================================================================
//
//
void print_db_file( const char *file ){
	
	printf("================================\n");
	printf("contents of %s\n",file);

	printf("theta folder            : %s \n", db.theta_folder );
	printf("theta file              : %s \n", db.theta_file );
	printf("logev file              : %s \n", db.logev_file );
	printf("theta folder            : %s \n", db.psi_folder );
	printf("psi file                : %s \n", db.psi_file );

	printf("Number of theta samples : %d \n", db.Ntheta );
	printf("Number of theta dim.    : %d \n", db.Dtheta );
	printf("Number of psi   samples : %d \n", db.Npsi );
	printf("Number of psi   dim.    : %d \n", db.Dpsi );
	
	printf("\n");
	printf("================================\n");

}






//=============================================================================
//
//
void loglike_posterior_theta_initialize( ) {




	read_db_file("db_theta.par");
	print_db_file("db_theta.par");
	
	// allocate global and local arrays
    double **theta    = (double **)malloc( sizeof(double *) * db.Ntheta );
    theta[0] = (double * )malloc( sizeof(double) * db.Ntheta * db.Dtheta);
 	for( int i=0; i<db.Ntheta; i++)
        theta[i] = ( *theta + db.Dtheta * i);

    double **psi    = (double **)malloc( sizeof(double *) * db.Npsi );
    psi[0] = (double * )malloc( sizeof(double) * db.Npsi * db.Dpsi);
 	for( int i=0; i<db.Npsi; i++)
        psi[i] = ( *psi + db.Dpsi * i);

	double *theta_logprior = (double *) malloc( db.Ntheta * sizeof(double) );
	
	denom          = (double *) malloc( db.Npsi   * sizeof(double) );




    char filename[BUFLEN];
    FILE *fp;

    // 1. read theta
    printf("\n1) Reading %d data from theta database.\n", db.Ntheta);

    snprintf(filename, BUFLEN, "%s/%s", db.theta_folder, db.theta_file );
    fp = fopen(filename, "r");
    if( fp == NULL ){
        printf("%s does not exist. Exit...\n", filename);
        exit(EXIT_FAILURE);
    }

    if( !enough_lines(fp, db.Ntheta) ){
            printf("\n\nError: Number of samples less than %d in file %s. Exit... \n\n", db.Ntheta, filename);
            exit(1);
    }

    for( int i=0; i<db.Ntheta; i++){
        for( int j=0; j<db.Dtheta; j++){
            fscanf(fp, "%lf", &theta[i][j]);
        }    double ignore;
        fscanf(fp, "%lf", &ignore);
        fscanf(fp, "%lf", &theta_logprior[i]);
    }
    fclose(fp);

    printf("\nSuccesfull reading data from theta database.\n\n");



    // 2. read log-evidence
    printf("\n2) Reading log-evidence.  \n");

	snprintf(filename, BUFLEN, "%s/%s", db.theta_folder, db.logev_file );
    fp = fopen(filename, "r");
    if (fp == NULL){
        printf("%s does not exist! exiting...\n", filename);
        exit(EXIT_FAILURE);
    }

    if (!enough_lines(fp, 1)){
    	printf("\n\n Error: Number of records less than %d in file %s. Exit... \n\n", 1, filename);
  		exit(EXIT_FAILURE);
    }

    fscanf(fp, "%lf", &logev_theta);
    fclose(fp);

    printf("\nSuccesfull reading of log-evidence: %.16f.\n\n", logev_theta);



    // 3a. read psi
    printf("\n3) Reading %d data from psi database. \n", db.Npsi);

    snprintf(filename, BUFLEN, "%s/%s", db.psi_folder, db.psi_file );
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s does not exist! exiting...\n", filename);
        exit( EXIT_FAILURE );
    }

    if( !enough_lines(fp, db.Npsi) ){
		printf("\n\n Error: Number of samples less than %d in file %s. Exit... \n\n", db.Npsi, filename);
 		exit(EXIT_FAILURE);
    }

    for( int i=0; i < db.Npsi; i++) {
        for( int j=0; j<db.Dpsi; j++)
            fscanf(fp, "%lf", &psi[i][j]);
        double ignore;
        fscanf(fp, "%lf", &ignore);
		fscanf(fp, "%lf", &ignore);
    }
    fclose(fp);
    printf("\nSuccesfull reading data from psi database.\n\n");



	// 3b. populate the priors
	priors = (Density **) malloc( db.Npsi*sizeof(Density*) );
	
	read_priors( "priors_theta.par", &priors[0], &Npr );
	reassign_prior( priors[0], Npr, psi[0] );
	
	for( int i=1; i<db.Ntheta; i++ ){
		new_prior_from_prior( &priors[i], priors[0], Npr );
		reassign_prior( priors[i], Npr, psi[i] );
	}

	




    // 4. compute denom
    printf("\n4) Computing the denominator....\n");

    double A[db.Ntheta];
    double B[db.Ntheta];  // exp(theta_logprior)

    for( int i=0; i<db.Npsi; i++)	denom[i] = 0;

    for( int i=0; i<db.Ntheta; i++)	B[i] = exp( theta_logprior[i] );

    for( int i = 0; i<db.Npsi ; i++){
        
		double sum=0;
        for( int j=0; j<db.Ntheta; j++){
            
            A[j] = prior_pdf( priors[i], Npr, theta[j]);
            sum +=  A[j] / B[j];
            
            //for( int k=0; k<db.Dpsi; k++)
            //    printf(" %lf  ", psi[j][k] );
            //printf("\n");
            //for( int k=0; k<db.Dtheta; k++)
            //    printf("    %lf  ", theta[j][k] );
            //printf("%lf - ", A[j] );
        }

        //printf(" \n");

        denom[i] = ( db.Npsi / db.Ntheta ) * sum;


    }

	printf("\nSuccesfull computation of the denominator.\n\n");

    //exit(1);



	free(theta[0]); free(theta);
	free(psi[0]); 	free(psi);
	free(theta_logprior);
}








//=============================================================================
//
//
void loglike_posterior_theta_finalize(){

	for(int i=0; i<db.Npsi; i++)
		delete_prior( priors[i], Npr );
	free(priors);


	free(denom);
}






//=============================================================================
//
//
double loglike_posterior_theta(double *theta, int n, void *output, int *info) {

	double out = 0;
	double sum = 0;

	for (int i = 0; i < db.Npsi; i++){
		double pr_hb  = exp( prior_log_pdf( priors[i], Npr, theta) );
		sum += pr_hb/denom[i];
	}


	if(sum==0)	return -1e12;

	double loglike_theta = loglike_(theta, n, output, info);
	if (loglike_theta == -1e12) {
		return loglike_theta;
	}

	out = loglike_theta - logev_theta + log(sum);

	if (isinf(out) || isnan(out)) out = -1e12;


	return out;

}




