/*
 *  engine_tmcmc.h
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _ENGINE_TMCMC_H_
#define _ENGINE_TMCMC_H_

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#ifdef __cplusplus
#include <mpi.h>
extern "C"
{
#endif
#include <torc.h>

#ifdef __cplusplus
}
#endif

#include <unistd.h>

#include "gsl_headers.h"

#define EXPERIMENTAL_RESULTS    0

/*** HELPER STRUCTS ***/
typedef struct data_s {
    int    Nth;        /* = PROBDIM*/
    int    MaxStages;    /* = MAXGENS*/
    int    PopSize;    /* = DATANUM*/

    double    *lowerbound;    /*[PROBDIM];*/
    double    *upperbound;    /*[PROBDIM];*/

    double *compositeprior_distr; /*[PROBDIM]*/

    double *prior_mu;
    double *prior_sigma;

    int auxil_size;
    double *auxil_data;

    int MinChainLength, MaxChainLength;

    double lb, ub;        /*generic lower and upper bound*/

    double    TolCOV;
    double    bbeta;
    long    seed;
    int    burn_in;

    struct optim_options {
        int    MaxIter;
        double    Tol;
        int    Display;
        double  Step;
    } options;

#if 1
    int    prior_type;     /* 0: uniform, 1: gaussian, 3: composite */
    int    load_from_file;
#endif

    int    icdump;
    int    ifdump;

    int    *Num;        /*[MAXGENS];*/
    int    LastNum;

    int use_proposal_cma;
    double  **init_mean;    /* [DATANUM][PROBDIM] */

    double  **local_cov;    /* [DATANUM][PROBDIM*PROBDIM] */
    int use_local_cov;
    double local_scale;

    int stealing;
    int restart;
} data_t;

typedef struct runinfo_s {
    int     Gen;
    double    *CoefVar;        /*[MAXGENS];*/
    double    *p;            /*[MAXGENS];        // cluster-wide*/
    int    *currentuniques;    /*[MAXGENS];*/
    double    *logselection;        /*[MAXGENS];*/
    double    *acceptance;        /*[MAXGENS];*/
    double    **SS;            /*[PROBDIM][PROBDIM];    // cluster-wide*/
    double    **meantheta;         /*[MAXGENS][PROBDIM]*/
} runinfo_t;

typedef struct {
    int idx;
    int nsel;
    double F;
} sort_t;

/*** DATABASES ***/
typedef struct cgdbp_s {
    double *point; /*[PROBDIM];*/
    double F;
	double prior;

    int counter;    /* not used (?)*/
    int nsel;    /* for selection of leaders only*/
    int queue;    /* for submission of leaders only*/
#if 1   // NN
    int surrogate;
    double error;
#endif
} cgdbp_t;

typedef struct cgdb_s {
    cgdbp_t *entry; /*[MAX_DB_ENTRIES];*/
    int entries;
    pthread_mutex_t m;
} cgdb_t;

typedef struct dbp_s {
    double *point; /*[PROBDIM];*/
    double F;
    int nG;
    double G[64];    /* maxG*/
    int surrogate;
} dbp_t;

typedef struct db_s {
    dbp_t *entry; /*[MAX_DB_ENTRIES];*/        /* */
    int entries;
    pthread_mutex_t m;
} db_t;

typedef struct resdbp_s {
    double *point;    /*[EXPERIMENTAL_RESULTS+1]; // +1 for the result (F)*/
    double F;
    int counter;    /* not used (?)*/
    int nsel;    /* for selection of leaders only*/
} resdbp_t;

typedef struct resdb_s {
    resdbp_t *entry; /*[MAX_DB_ENTRIES];*/
    int entries;
    pthread_mutex_t m;
} resdb_t;
/*** END HELPER STRUCTS ***/

/*** DATABASE INSTANCES ***/
extern data_t data;
extern runinfo_t runinfo;
extern cgdb_t curgen_db;
extern db_t full_db;
extern resdb_t curres_db;

void update_full_db(double point[], double F, double *G, int n, int surrogate);
void init_full_db();

void update_curgen_db(double point[], double F, double prior);
void init_curgen_db();

void update_curres_db(double point[], double F);
void init_curres_db();
void print_full_db();
void print_curgen_db();
void dump_curgen_db(int Gen);
void dump_curres_db(int Gen);
void dump_full_db(int Gen);
void display_curgen_db(int Gen);
int load_curgen_db(int Gen);

/*** UTILS ***/
double compute_sum(double *x, int n);
double compute_mean(double *x, int n);
double compute_std(double *x, int n, double mean);
double compute_min(double *x, int n);
int compute_min_idx_i(int *v, int n);
double compute_max(double *x, int n);
void print_matrix(char *name, double *x, int n);
void print_matrix_i(char *name, int *x, int n);
void print_matrix_2d(char *name, double **x, int n1, int n2);

/*** RNG ***/
void gsl_rand_init(int seed);
double normalrand(double mu, double var);
double uniformrand(double a, double b);
void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[]);
void shuffle(int *perm, int N);
int mvnrnd(double *mean, double *var, double *res, int n);
double mvnpdf(int n, double *xv, double *mv, double *vm);
double logmvnpdf(int n, double *xv, double *mv, double *vm);

double truncated_normal_pdf (double x, double mu, double sigma, double a, double b);
double truncated_normal_rand (double mu, double sigma, double a, double b);
double truncated_lognormal_pdf (double x, double mu, double sigma, double a, double b);

/*** STATISTICS ***/
void calculate_statistics(double flc[], int n, int nselections, int gen, unsigned int sel[]);

/*** PROBLEM FUNCTIONS ***/
double likelihood(double *x, int N);
double posterior(double *theta, int n, double LH);
double logpriorpdf(double *theta, int n);

/*** AUX ***/
void inc_nfc();
void get_nfc_task(int *);
int get_nfc();
void reset_nfc_task();
void reset_nfc();
int get_tfc();

/*** POSDEF ***/
void compute_mat_product_vect(double *mat/*2D*/, double vect[], double res_vect[], double coef, int PROBDIM);
double compute_dot_product(double row_vector[], double vector[], int PROBDIM);
int inv_matrix(double coef, double *current_hessian/*2D*/, double *inv_current_hessian/*2D*/, int PROBDIM);

#endif
