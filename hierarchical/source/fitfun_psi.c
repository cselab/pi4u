#include "fitfun.h"

#include "gsl_headers.h"


#define DATABASE  "../../data/theta"
#define THETAFILE1 "theta_"
#define THETAFILE2 ".txt"
#define LOGEVFILE1 "evidence_"
#define LOGEVFILE2 ".txt"

#define BUFLEN 1024

#define NREC 1000

#define N_IND 5

static int ind_list[N_IND] = {1, 2, 3, 4, 5};



typedef struct ffdata_s {
    double *x;
    double loglik;
    double logprior;
} ffdata_t;

static ffdata_t *ffdata[N_IND];
static double logEv[N_IND];







double loglike_psi(double *psi, int n);


#include "engine_tmcmc.h"
extern data_t data;

void fitfun_initialize() {

    int nn = data.Nth;

    printf("\nReading %d data from theta database for individuals: \n", NREC);
    for (int i = 0; i < N_IND; i++)
        printf("%d  -  ", ind_list[i]);
    printf("\n");

    int n = nn/2;
    for (int pp = 1; pp <= N_IND; pp++) {
        int p = pp-1;
        char filename[BUFLEN];
        FILE *fp;


        snprintf(filename, BUFLEN, "%s/%s%03d%s", DATABASE, THETAFILE1,  ind_list[p], THETAFILE2);
        fp = fopen(filename, "r");
        if (fp == NULL) {
            printf("\n %s  does not exist. Exiting...\n", filename);
            exit(1);
        }

        // count the number of lines in file and check if less than NREC
        char ch;
        int lines = 0;

        while (!feof(fp)) {
            ch = fgetc(fp);
            if (ch == '\n')
                lines++;
        }
        rewind(fp);

        if (lines < NREC) {
            printf("\n\n Error: Number of samples less than %d in file %s. "
                    "Exit... \n\n", NREC, filename);
            exit(1);
        }

        ffdata[p] = (ffdata_t*)malloc(NREC*sizeof(ffdata_t));
        for (int i = 0; i < NREC; i++)
            ffdata[p][i].x = (double*)malloc(n*sizeof(double));

        for (int i = 0; i < NREC; i++) {
            for (int j = 0; j < n; j++)
                fscanf(fp, "%lf", &(ffdata[p][i].x[j]));

            fscanf(fp, "%lf", &ffdata[p][i].loglik);
            fscanf(fp, "%lf", &ffdata[p][i].logprior);

            // ffdata[p][i] = exp(ffdata[p][i].logprior);
        }
        fclose(fp);

        snprintf(filename, BUFLEN, "%s/%s%03d%s", DATABASE, LOGEVFILE1,  ind_list[p], LOGEVFILE2);
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


void fitfun_finalize() {
    for (int pp = 1; pp <= N_IND; pp++) {
        int p = pp-1;
        if (ffdata[p] != NULL) {
            for (int i = 0; i < NREC; i++) free(ffdata[p][i].x);
            free(ffdata[p]);
        }
    }
}




double fitfun(double *x, int n, void *output, int *info) {

    double res;
    res = loglike_psi(x, n);

    int BUFLEN_TMP=64;
    char msg[BUFLEN], tmp[BUFLEN_TMP];
    snprintf(msg, BUFLEN, "task: ");

    for (int i = 0; i < n; i++) {
        snprintf(tmp, BUFLEN_TMP, "%lf ", x[i]);
        snprintf(msg, BUFLEN, "%s %s", msg, tmp);
    }

    snprintf(tmp, BUFLEN_TMP, "%lf\n", res);
    snprintf(msg, BUFLEN, "%s %s", msg, tmp);

    return res;
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



inline double trnpdf(double x, double m, double s, double l, double u) {
    if (s > 0)
        return gsl_ran_gaussian_pdf(fabs(x-m), s) /
            (1.0 - gsl_cdf_gaussian_P(fabs(m-l), s)
                 - gsl_cdf_gaussian_P(fabs(m-u), s));
    else
        return NAN;
}





inline double log_lognpdf(double x, double m, double s) {
    if (s > 0)
      return -0.5*log(2*M_PI) - log(x) - log(s) -0.5*pow((log(x)-m)/s,2);
    else
      return NAN;
}



double log_priorHB(double *x, double *psi, int n) {

    double res = 0.0;

    for (int i = 0; i < n/2; i++) {
        res += log_lognpdf(x[i], psi[2*i], psi[2*i+1]);
    }

    return res;
}




inline double loglike_psi(double *psi, int n) {
    double out = 0;
            
	//for(int k=0; k<n; k++) printf("%f \t ",psi[k]);
    //printf("\n");

    for (int pp = 1; pp <= N_IND; pp++) {
        int p = pp-1;

        if (ffdata[p] == NULL) {
            printf("\n Pointer in %d file  is NULL. Exiting...\n", p);
            exit(1);
        }

        double sum = 0;
        for (int i = 0; i < NREC; i++) {

            double log_pr_hb = log_priorHB(ffdata[p][i].x, psi, n);
            double log_pri = ffdata[p][i].logprior;

            // for(int k=0; k<n/2; k++) printf("%f \t ",ffdata[p][i].x[k]);
            // printf("\n");
            // printf("%f \n",log_pr_hb);
            // printf("%f \n",log_pri);
            // printf("%f \n",log_pr_hb-log_pri);
            // exit(0);

            sum += exp( log_pr_hb-log_pri ) ;
        }

        if (sum == 0 || isnan(sum) || isinf(sum)) {
            out = -1e12;
            return out;
        }

        out = out + logEv[p] - log(NREC) + log(sum);
    }

    if (isinf(out) || isnan(out)) out = -1e12;

    return out;
}
