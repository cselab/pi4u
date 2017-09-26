/*
 *  engine_tmcmc.h
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
//#include <mpi.h>
//#include <torc.h>

#include "gsl_headers.h"

typedef struct data_s {
	int	Nth;		// = PROBDIM
	int	MaxStages;	// = MAXGENS
	int	PopSize;	// = DATANUM

	double	*lowerbound;	//[PROBDIM];
	double	*upperbound;	//[PROBDIM];

	double lb, ub;		//generic lower and upper bound

	double	TolCOV;
	double	bbeta;
	int	seed;

	struct optim_options {
		int	MaxIter;
		double	Tol;
		int	Display;
	} options;

	int	iplot;
	int	idump;

	int	*Num;		//[MAXGENS];
	int	LastNum;

} data_t;

typedef struct runinfo_s {
	int 	Gen;
	double	*CoefVar;		//[MAXGENS];
	double	*p;			//[MAXGENS];		// cluster-wide
	int	*currentuniques;	//[MAXGENS];
	double	*logselection;		//[MAXGENS];
	double	*acceptance;		//[MAXGENS];
	double	**SS;			//[PROBDIM][PROBDIM];	// cluster-wide
	double	**meantheta; 		//[MAXGENS][PROBDIM]
} runinfo_t;


/*** DATABASES ***/
typedef struct cgdbp_s {
	double *point; //[PROBDIM];
	double F;
	int counter;	// not used (?)
	int nsel;	// for selection of leaders only
	int queue;	// for submission of leaders only
} cgdbp_t;

typedef struct cgdb_s {
	cgdbp_t *entry; //[MAX_DB_ENTRIES];
	int entries;
#if 0
	pthread_mutex_t m;
#endif
} cgdb_t;


typedef struct tmcmc_data_s {
	data_t data;
	runinfo_t runinfo;
	cgdb_t curgen_db;
} tmcmc_data_t;

extern tmcmc_data_t tmcmc_data;
//extern data_t data;
//extern runinfo_t runinfo;
//extern cgdb_t curgen_db;

void update_curgen_db(cgdb_t *curgen_db, double point[], double F, int dim);
void init_curgen_db(cgdb_t *curgen_db, int popsize);
void dump_curgen_db(cgdb_t *curgen_db, int Gen, int dim, int tinfo[4]);

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

/*** STATISTICS ***/
void calculate_statistics(tmcmc_data_t *tmcmc_data, double flc[], int n, int nselections, int gen, unsigned int sel[]);

/*** PROBLEM FUNCTIONS ***/
double likelihood(double *x, int N);
double posterior(double *theta, int n, double LH);
double priorpdf(double *theta, int n, double *lowerbound, double *upperbound);

/*** AUX ***/
void inc_nfc(void);
void get_nfc_task(int *);
int get_nfc(void);
void reset_nfc_task(void);
void reset_nfc(void);
int get_tfc(void);


int tmcmc_initialize(char *fitfun_name);
void tmcmc_finalize(void);
void tmcmc(double *res, int tmcmc_info[4], int Nth, int MaxStages, int PopSize, double *lb, double *ub);
