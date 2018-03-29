#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gsl_headers.h"
#include "fitfun_argon.c"

#define DATABASE   "/cluster/home02/mavt/kulina/UQ/engines_argon/hierarchical/cpp/database/p_free"
#define THETAFILE  "final.txt"
#define NTHETA     4000
#define DIMTHETA   4
#define LOGEVFILE  "fitness.txt"

#define DATASETNUM 1  // number of the dataset for which we are predicting

#define PSIFILE    "psi_lv_unif_a_l.txt"
#define NPSI       60000

#define BUFLEN 1024

static double theta[NTHETA][DIMTHETA];
static double psi[NPSI][DIMTHETA*2];
static double logev_theta;
static double theta_logprior[NTHETA];
static double denom[NPSI];



// count the number of lines in file fp  and check if less than N
int enough_lines(FILE *fp, int N) {
    char ch;
    int lines = 0;

    while (!feof(fp)) {
            ch = fgetc(fp);
            if (ch == '\n')
                lines++;
    }

    rewind(fp);

    return lines >= N;
}



inline double unifpdf(double x, double l, double u) {
    if ((l <= x) && (x <= u))
        return 1.0/(u-l);
    return 0;
}



inline double unifpdf2(double x, double l, double len) {
    if ((l <= x) && (x <= l+len))
        return 1.0/len;
    return 0;
}



inline double lognpdf(double x, double m, double s) {
    if (s > 0)
        return gsl_ran_lognormal_pdf(x, m, s);
    else
        return NAN;
}


inline double trnpdf(double x, double m, double s, double l, double u) {
    if (s > 0)
        return gsl_ran_gaussian_pdf(fabs(x-m), s) /
            (1.0 - gsl_cdf_gaussian_P(fabs(m-l), s)
                 - gsl_cdf_gaussian_P(fabs(m-u), s));
    else
        return NAN;
}



inline double priorHB(double *x, double *psi, int n) {
    double res = 1;

    /* for the truncated gaussian */
    double l[4] = {0.05, 3., 6, 1e-6};
    double u[4] = {10.0, 4., 15, 1};

    for (int i = 0; i < n/2; i++) {
        /* res *= unifpdf(x[i], psi[2*i], psi[2*i+1]); */
        res *= unifpdf2(x[i], psi[2*i], psi[2*i+1]);
        /* res *= lognpdf(x[i]-l[i], psi[2*i], psi[2*i+1]); */
        /* res *= trnpdf(x[i], psi[2*i], psi[2*i+1], l[i], u[i]); */
    }

    return res;
}







void fitfun_init(int n) {
    char filename[BUFLEN];
    FILE *fp;

    // 1. read theta
    printf("\n1) Reading %d data from theta database for individual: %d  \n",
            NTHETA, DATASETNUM);

    snprintf(filename, BUFLEN, "%s/tmcmc_dir%d/%s", DATABASE, DATASETNUM,
            THETAFILE);
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s does not exist! exiting...\n", filename);
        exit(1);
    }

    if (!enough_lines(fp, NTHETA)) {
            printf("\n\n Error: Number of samples less than %d in file %s. "
                    "Exit... \n\n", NTHETA, filename);
            exit(1);
    }

    for (int i = 0; i < NTHETA; i++) {
        for (int j = 0; j < n; j++) {
            fscanf(fp, "%lf", &theta[i][j]);
        }
        double ignore;
        fscanf(fp, "%lf", &ignore);
        fscanf(fp, "%lf", &theta_logprior[i]);
    }
    fclose(fp);

    printf("\nSuccesfull reading data from theta database.\n\n");

    // 2. read fitness
    printf("\n2) Reading log-evidence.  \n");

    snprintf(filename, BUFLEN, "%s/tmcmc_dir%d/%s", DATABASE, DATASETNUM,
            LOGEVFILE);
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s does not exist! exiting...\n", filename);
        exit(1);
    }

    if (!enough_lines(fp, 1)) {
            printf("\n\n Error: Number of records less than %d in file %s. "
                    "Exit... \n\n", 1, filename);
            exit(1);
    }

    fscanf(fp, "%lf", &logev_theta);
    fclose(fp);

    printf("\nSuccesfull reading of log-evidence: %.16f.\n\n", logev_theta);

    // 3. read psi
    printf("\n3) Reading %d data from psi database. \n", NPSI);

    snprintf(filename, BUFLEN, "%s/%s", DATABASE, PSIFILE);
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("%s does not exist! exiting...\n", filename);
        exit(1);
    }

    if (!enough_lines(fp, NPSI)) {
            printf("\n\n Error: Number of samples less than %d in file %s. "
                    "Exit... \n\n", NTHETA, filename);
            exit(1);
    }


    for (int i = 0; i < NPSI; i++) {
        for (int j = 0; j < 2*n; j++)
            fscanf(fp, "%lf", &psi[i][j]);
        double ignore;
        fscanf(fp, "%lf", &ignore);
    }
    fclose(fp);
    printf("\nSuccesfull reading data from psi database.\n\n");

    // 4. compute denom
    printf("\n4) Computing the denominator....\n");

    static double B[NTHETA];  // exp(theta_logprior)
    static double A[NTHETA];

    double t0 = omp_get_wtime();
    for (int i = 0; i < NPSI; i++)
        denom[i] = 0;

    for (int i = 0; i < NTHETA; i++)
        B[i] = exp(theta_logprior[i]);

#pragma omp parallel for
    for (int i = 0; i < NPSI; i++) {
        double sum;
        sum = 0;
        for (int r = 0; r < NTHETA; r++) {
            A[r] = priorHB(theta[r], psi[i], n);
            sum += (A[r]/B[r]);
        }
        denom[i] = (NPSI/NTHETA)*sum;
    }
    double t1 = omp_get_wtime();
    printf("\nSuccesfull computation of the denominator in %f seconds.\n\n",
            t1-t0);
}








void fitfun_cleanup() {
}











double posterior_theta(double *theta, int n, void *output, int *info) {
    double out = 0;
    double sum = 0;

    for (int i = 0; i < NPSI; i++){
        double pr_hb  = priorHB(theta, psi[i], n);
        sum += pr_hb/denom[i];
    }

    if (sum == 0)	return = -1e12;



    double loglike_theta = fitfun(theta, n, output, info);
	if (loglike_theta == -1e12) {
		return loglike_theta;
	}

    out = loglike_theta - logev_theta + log(sum);

	printf("posterior_theta: ( ");
	for (int i = 0; i < n; ++i)
		printf("%.16f ", theta[i]);
	printf(") = {%.16f, %.16f, %.16f}\n", loglike_theta, -logev_theta, log(sum));

    if (isinf(out) || isnan(out)) out = -1e12;

    return out;
}




double fitfun(double *x, int n, void *output, int *info) {
    return posterior_theta(x, n, output, info);
}
