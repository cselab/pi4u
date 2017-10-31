#include "fitfun.h"

#include "gsl_headers.h"



#define DATABASE  "./data/theta"
#define THETAFILE1 "theta_"
#define THETAFILE2 ".txt"
#define LOGEVFILE1 "evidence_"
#define LOGEVFILE2 ".txt"

#define BUFLEN 1024

#define NREC 10000

#define NPAT 5

static int pat_list[NPAT] = { 1, 2, 3, 4, 5 };


typedef struct ffdata_s {
    double *x;
    double loglik;
    double prior;
} ffdata_t;

static ffdata_t *ffdata[NPAT];
static double logEv[NPAT];







double loglike_psi(double *psi, int n);



void fitfun_init( data_t data ) {

    int nn = data.Nth;

    printf("\nReading %d data from theta database for individuals: \n", NREC);
    for (int i = 0; i < NPAT; i++)
        printf("%d  -  ", pat_list[i]);
    printf("");

    int n = nn/2;
    for (int pp = 1; pp <= NPAT; pp++) {
        int p = pp-1;
        char filename[BUFLEN];
        FILE *fp;


        snprintf(filename, BUFLEN, "%s/%s%03d%s", DATABASE, THETAFILE1,  pat_list[p], THETAFILE2);
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
            fscanf(fp, "%lf", &ffdata[p][i].prior);

            ffdata[p][i].prior = exp(ffdata[p][i].prior);
        }
        fclose(fp);

        snprintf(filename, BUFLEN, "%s/%s%03d%s", DATABASE, LOGEVFILE1,  pat_list[p], LOGEVFILE2);
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


void fitfun_cleanup() {
    for (int pp = 1; pp <= NPAT; pp++) {
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


double priorHB(double *x, double *psi, int n) {
    double res = 1.0;

    for (int i = 0; i < n/2; i++) {
        res *= lognpdf(x[i], psi[2*i], psi[2*i+1]);
    }


    return res;

}


inline double loglike_psi(double *psi, int n) {
    double out = 0;

    for (int pp = 1; pp <= NPAT; pp++) {
        int p = pp-1;

        if (ffdata[p] == NULL) {
            printf("\n Pointer in %d file  is NULL. Exiting...\n", p);
            exit(1);
        }

        double sum = 0;
        for (int i = 0; i < NREC; i++) {
            double pr_hb = priorHB(ffdata[p][i].x, psi, n);
            double pri = ffdata[p][i].prior;
            sum += pr_hb/pri;
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
