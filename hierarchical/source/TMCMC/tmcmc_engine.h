/*
 *  tmcmc_engine.h
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _ENGINE_TMCMC_H_
#define _ENGINE_TMCMC_H_




#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>



#if defined(_USE_TORC_)
	#include <mpi.h>
	#ifdef __cplusplus
		extern "C"
		{
	#endif

	#include <torc.h>

	#ifdef __cplusplus
		}
	#endif

#else
	#include <pthread.h>
#endif



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

    int    prior_type;     /* 0: uniform, 1: gaussian, 3: composite */
    int    load_from_file;

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
    int surrogate;
    double error;
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






void read_data();
void data_init();
void save_runinfo();
int load_runinfo();
void check_for_exit();
void torc_update_full_db_task(double point[], double *pF, int *psurrogate);
void torc_update_full_db_task_p5(double point[], double *pF, double *G, int *pn, int *psurrogate);
void torc_update_full_db(double point[], double F, double *G, int n, int surrogate);
void torc_update_curgen_db_task(double point[], double *pF, double *pprior);
void torc_update_curgen_db(double point[], double F, double prior);
void torc_update_curres_db_task(double point[EXPERIMENTAL_RESULTS], double *pF);
void torc_update_curres_db(double point[EXPERIMENTAL_RESULTS], double F);

/*** TASK MANAGEMENT ***/
void taskfun(double /*const*/ *x, int *pN, double *res, int winfo[4]);
double F(double *TP, int *pn);    /* for PNDL */
void evaluate_F(double point[], double *Fval, int worker_id, int gen_id, int chain_id, int step_id, int ntasks);
void initchaintask(double in_tparam[], int *pdim, double *out_tparam, int winfo[4]);
static int in_rect(double *v1, double *v2, double *diam, double sc, int D);
void precompute_chain_covariances(const cgdbp_t* leader,double** init_mean, double** chain_cov, int newchains);
int compute_candidate(double candidate[], double chain_mean[], double var);
int compute_candidate_cov(double candidate[], double chain_mean[], double chain_cov[]);
void chaintask(double in_tparam[], int *pdim, int *pnsteps, double *out_tparam, int winfo[4],double *init_mean, double *chain_cov);
int compar_desc(const void* p1, const void* p2);
int prepare_newgen(int nchains, cgdbp_t *leaders);







#endif
