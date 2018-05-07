#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>

#include "loglike_psi.h"
#include "../priors/priors.h"

#define BUFLEN 1024

#define BLANKS " \t\n"


typedef struct theta_database_s{
    char folder[128];
    char theta_prefix[128];
    char theta_suffix[128];
    char logev_prefix[128];
    char logev_suffix[128];
    int  Nrec;
    
    int Nind;
    int *ind_list;
}theta_database;


typedef struct theta_data_s{
    double *x;
    double loglik;
    double logprior;
} theta_data_t;



static theta_database  theta_db;

static theta_data_t **theta_data;
static double *logEv;

static Density *prior_theta=NULL;
static int Npr=0;









//=============================================================================
//
//
int count_lines(FILE *fp){

    char ch;
    int lines = 0;
    
    rewind(fp);
    
    while (!feof(fp)) {
        ch = fgetc(fp);
        if (ch == '\n')
            lines++;
    }
    
    rewind(fp);

    return lines;

}




//=============================================================================
//
//
void read_theta_file( const char *file  ){

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
            
            if( strcmp(pch,"folder")==0 ){
                pch = strtok (NULL, BLANKS );
                strcpy( theta_db.folder, pch );
                break;
            }

            if( strcmp(pch,"theta_prefix")==0 ){
                pch = strtok (NULL, BLANKS );
                strcpy( theta_db.theta_prefix, pch );
                break;
            }

            if( strcmp(pch,"theta_suffix")==0 ){
                pch = strtok (NULL, BLANKS );
                strcpy( theta_db.theta_suffix, pch );
                break;
            }

            if( strcmp(pch,"logev_prefix")==0 ){
                pch = strtok (NULL, BLANKS );
                strcpy( theta_db.logev_prefix, pch );
                break;
            }

            if( strcmp(pch,"logev_suffix")==0 ){
                pch = strtok (NULL, BLANKS );
                strcpy( theta_db.logev_suffix, pch );
                break;
            }

            if( strcmp(pch,"Nrec")==0 ){
                pch = strtok (NULL, BLANKS );
                theta_db.Nrec = atoi( pch );
                break;
            }

            if( strcmp(pch,"Nind")==0 ){
                pch = strtok (NULL, BLANKS );
                theta_db.Nind = atoi( pch );
                theta_db.ind_list = (int *) malloc( theta_db.Nind*sizeof(int) );
                break;
            }

            if( strcmp(pch,"ind_list")==0 ){
                if( theta_db.Nind==0 || theta_db.ind_list==NULL  ){
                    printf("Define 'Nind' before 'ind_list' in %s. Exit...",file);
                    exit(EXIT_FAILURE);
                }
                
                for(int i=0; i<theta_db.Nind; i++){ //XXX: more checks here
                    pch = strtok (NULL, BLANKS );
                    theta_db.ind_list[i] = atoi( pch );
                }
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
void print_theta_file( const char *file ){
    
    printf("================================\n");
    printf("contents of %s\n",file);

    printf("data base folder        : %s \n",theta_db.folder);
    printf("theta files             : %s*%s \n",theta_db.theta_prefix, theta_db.theta_suffix);
    printf("logev files             : %s*%s \n",theta_db.logev_prefix, theta_db.logev_suffix);
    printf("Number of samples       : %d \n", theta_db.Nrec );
    printf("Number of individuals   : %d \n", theta_db.Nind );
    printf("Individuals             : ");
    for(int i=0; i<theta_db.Nind; i++)
        printf("%d -  ",theta_db.ind_list[i]);  
    
    printf("\n================================\n");

}




//=============================================================================
//
//
void loglike_psi_initialize( int nn ){

    read_priors( "priors_theta.par", &prior_theta, &Npr );

    read_theta_file( "db_psi.par" );
    print_theta_file( "db_psi.par" );


    printf("\nReading %d data from theta database for individuals: \n", theta_db.Nrec);
    for (int i = 0; i < theta_db.Nind; i++)
        printf("%d  -  ", theta_db.ind_list[i]);
    printf("\n");
    
    theta_data = (theta_data_t**) malloc( theta_db.Nind*sizeof(theta_data_t*) );
    logEv = (double *) malloc( theta_db.Nind*sizeof(double) );


    for( int p=0; p < theta_db.Nind; p++ ){
        
        char filename[BUFLEN];
        FILE *fp;


        snprintf(filename, BUFLEN, "%s/%s%03d%s", theta_db.folder, theta_db.theta_prefix,  theta_db.ind_list[p], theta_db.theta_suffix);
        fp = fopen(filename, "r");
        if (fp == NULL) {
            printf("\n %s  does not exist. Exiting...\n", filename);
            exit(1);
        }

        if ( count_lines(fp) < theta_db.Nrec){
            printf("\n\n Error: Number of samples less than %d in file %s.Exit... \n\n", theta_db.Nrec, filename);
            exit(EXIT_FAILURE);
        }

        theta_data[p] = (theta_data_t*)malloc( theta_db.Nrec*sizeof(theta_data_t) );
        for (int i = 0; i < theta_db.Nrec; i++)
            theta_data[p][i].x = (double*)malloc( Npr*sizeof(double));

        for (int i = 0; i < theta_db.Nrec; i++) {
            for (int j = 0; j < Npr; j++)
                fscanf(fp, "%lf", &(theta_data[p][i].x[j]));

            fscanf(fp, "%lf", &theta_data[p][i].loglik);
            fscanf(fp, "%lf", &theta_data[p][i].logprior);
        }
        
        fclose(fp);

        snprintf(filename, BUFLEN, "%s/%s%03d%s", theta_db.folder, theta_db.logev_prefix,  theta_db.ind_list[p], theta_db.logev_suffix);
        fp = fopen(filename, "r");
        if (fp == NULL) {
            printf("\n %s  does not exist. Exiting...\n", filename);
            exit(1);
        }

        fscanf(fp, "%lf", &logEv[p]);
        fclose(fp);
    }

    printf("\nSuccesfull reading data from theta database.\n\n");
}











//=============================================================================
//
//
void loglike_psi_finalize() {

    printf("\nFinilizing psi loglike...\n");fflush(NULL);

    for( int p=0 ; p < theta_db.Nind;  p++ ){
        if (theta_data[p] != NULL){
            for (int i = 0; i < theta_db.Nrec; i++) free(theta_data[p][i].x);
            free(theta_data[p]);
        }
    }

    delete_prior( prior_theta, Npr );
}








//=============================================================================
//
//
double loglike_psi(double *psi, int n, void *output, int *info) {

    double out = 0;

    //printf("Nind %d \n",theta_db.Nind);
    //printf("psi: ");
    //for(int i=0; i<n; i++)
    //  printf("%lf ",psi[i]);
    //printf("\n");
    //printf("---------------\n");


    Density *prior_tmp=NULL;
    new_prior_from_prior(&prior_tmp, prior_theta, Npr);
    reassign_prior( prior_tmp, Npr, psi );
    
    //print_priors(prior_theta, Npr);

    for( int p=0; p<theta_db.Nind; p++){


        if( theta_data[p]==NULL ){
            printf("\n Pointer in %d file  is NULL. Exit...\n", p);
            exit(EXIT_FAILURE);
        }

        double sum = 0;
        for( int i = 0; i < theta_db.Nrec; i++ ){
    
            double log_pr_hb = prior_log_pdf( prior_tmp, Npr, theta_data[p][i].x );

            double log_pri = theta_data[p][i].logprior;

            //printf("      --->  ");
            //for(int j=0; j<Npr; j++)
            //  printf(" %lf ", theta_data[p][i].x[j] );
            //printf("  --- %lf ", log_pr_hb);
            //printf("\n");

            sum += exp( log_pr_hb-log_pri ) ;
        }


        //printf("%d -> %lf \n",p,sum);

        if( sum == 0 || isnan(sum) || isinf(sum) ){
            out = -1e12;
            return out;
        }

        out = out + logEv[p] - log(theta_db.Nrec) + log(sum);
    }

    if (isinf(out) || isnan(out)) out = -1e12;


    delete_prior(prior_tmp,Npr);


    //XXX: what is this? delete?
    int BUFLEN_TMP=64;
    char msg[BUFLEN], tmp[BUFLEN_TMP];
    snprintf(msg, BUFLEN, "task: ");

    for (int i = 0; i < n; i++) {
        snprintf(tmp, BUFLEN_TMP, "%lf ", psi[i]);
        snprintf(msg, BUFLEN, "%s %s", msg, tmp);
    }

    snprintf(tmp, BUFLEN_TMP, "%lf\n", out);
    snprintf(msg, BUFLEN, "%s %s", msg, tmp);
    //**********************************************


    return out;
}

















