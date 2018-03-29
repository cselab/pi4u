#include "loglike_psi.h"
#include "../TMCMC/engine_tmcmc.h"

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


extern data_t data;



void loglike_psi_initialize() {

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





void loglike_psi_finalize() {
    for (int pp = 1; pp <= N_IND; pp++) {
        int p = pp-1;
        if (ffdata[p] != NULL) {
            for (int i = 0; i < NREC; i++) free(ffdata[p][i].x);
            free(ffdata[p]);
        }
    }
}






double loglike_psi(double *psi, int n, void *output, int *info) {

	double out = 0;

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

			sum += exp( log_pr_hb-log_pri ) ;
		}

		if (sum == 0 || isnan(sum) || isinf(sum)) {
			out = -1e12;
			return out;
		}

		out = out + logEv[p] - log(NREC) + log(sum);
	}

	if (isinf(out) || isnan(out)) out = -1e12;


	int BUFLEN_TMP=64;
    char msg[BUFLEN], tmp[BUFLEN_TMP];
    snprintf(msg, BUFLEN, "task: ");

    for (int i = 0; i < n; i++) {
        snprintf(tmp, BUFLEN_TMP, "%lf ", psi[i]);
        snprintf(msg, BUFLEN, "%s %s", msg, tmp);
    }

    snprintf(tmp, BUFLEN_TMP, "%lf\n", out);
    snprintf(msg, BUFLEN, "%s %s", msg, tmp);

    return out;
}
