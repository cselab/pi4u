#include <fitfun.h>
#include "loglike_posterior_theta.h"
#include "priors.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>




#define BUFLEN 1024



#define THETADATABASE  "."
#define THETAFILE "theta.txt"
#define LOGEVFILE "evidence.txt"
#define DIMTHETA   4


#define PSIDATABASE  "."
#define PSIFILE "psi.txt"
#define DIMPSI	8


#define NTHETA     1000
#define NPSI       1000




static double theta[NTHETA][DIMTHETA];
static double psi[NPSI][DIMPSI];
static double logev_theta;
static double theta_logprior[NTHETA];
static double denom[NPSI];






// count the number of lines in file fp  and check if less than N
static int enough_lines(FILE *fp, int N) {
    char ch;
    int lines = 0;

    while (!feof(fp)){
            ch = fgetc(fp);
            if (ch == '\n')
                lines++;
    }
    rewind(fp);

    return lines >= N;
}








void loglike_posterior_theta_initialize( ) {

    char filename[BUFLEN];
    FILE *fp;

    // 1. read theta
    printf("\n1) Reading %d data from theta database.\n", NTHETA);

    snprintf(filename, BUFLEN, "%s/%s", THETADATABASE,THETAFILE);
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s does not exist! exiting...\n", filename);
        exit(1);
    }

    if (!enough_lines(fp, NTHETA)){
            printf("\n\n Error: Number of samples less than %d in file %s. Exit... \n\n", NTHETA, filename);
            exit(1);
    }

    for (int i = 0; i < NTHETA; i++){
        for (int j = 0; j < DIMTHETA; j++){
            fscanf(fp, "%lf", &theta[i][j]);
        }
        double ignore;
        fscanf(fp, "%lf", &ignore);
        fscanf(fp, "%lf", &theta_logprior[i]);
    }
    fclose(fp);

    printf("\nSuccesfull reading data from theta database.\n\n");



    // 2. read log-evidence
    printf("\n2) Reading log-evidence.  \n");

	snprintf(filename, BUFLEN, "%s/%s", THETADATABASE, LOGEVFILE);
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s does not exist! exiting...\n", filename);
        exit(1);
    }

    if (!enough_lines(fp, 1)){
            printf("\n\n Error: Number of records less than %d in file %s. Exit... \n\n", 1, filename);
            exit(1);
    }

    fscanf(fp, "%lf", &logev_theta);
    fclose(fp);

    printf("\nSuccesfull reading of log-evidence: %.16f.\n\n", logev_theta);



    // 3. read psi
    printf("\n3) Reading %d data from psi database. \n", NPSI);

    snprintf(filename, BUFLEN, "%s/%s", PSIDATABASE, PSIFILE);
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s does not exist! exiting...\n", filename);
        exit(1);
    }

    if (!enough_lines(fp, NPSI)) {
            printf("\n\n Error: Number of samples less than %d in file %s. Exit... \n\n", NPSI, filename);
            exit(1);
    }

    for (int i = 0; i < NPSI; i++) {
        for (int j = 0; j < DIMPSI; j++)
            fscanf(fp, "%lf", &psi[i][j]);
        double ignore;
        fscanf(fp, "%lf", &ignore);
		fscanf(fp, "%lf", &ignore);
    }
    fclose(fp);
    printf("\nSuccesfull reading data from psi database.\n\n");

    // 4. compute denom
    printf("\n4) Computing the denominator....\n");

    static double B[NTHETA];  // exp(theta_logprior)
    static double A[NTHETA];

    // double t0 = omp_get_wtime();
    for (int i = 0; i < NPSI; i++)	denom[i] = 0;

    for (int i = 0; i < NTHETA; i++) B[i] = exp(theta_logprior[i]);

// // #pragma omp parallel for
    for (int i = 0; i < NPSI; i++){
        double sum;
        sum = 0;
        for (int r = 0; r < NTHETA; r++){
            A[r] = log_priorHB(theta[r], psi[i], DIMTHETA) ;
            sum +=  A[r] / B[r];
        }
        denom[i] = ( NPSI / NTHETA ) * sum;

    }
    // double t1 = omp_get_wtime();
    // printf("\nSuccesfull computation of the denominator in %f seconds.\n\n", t1-t0);
	printf("\nSuccesfull computation of the denominator.\n\n");
}








void loglike_posterior_theta_finalize(){
}






double loglike_posterior_theta(double *theta, int n, void *output, int *info) {
	double out = 0;
	double sum = 0;

	for (int i = 0; i < NPSI; i++){
		double pr_hb  = exp( log_priorHB(theta, psi[i], DIMTHETA) );
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
