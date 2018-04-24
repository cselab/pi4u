/*
 *  tmcmc_aux.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <math.h>
#include <stdio.h>
#include <time.h>


#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>




#include "tmcmc_aux.h"
#include "tmcmc_engine.h"








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
	int torc_node_id() { return 0; }
	int torc_num_nodes() { return 1; }

	#if defined(_USE_OPENMP_)
		#include <omp.h>
		int torc_i_worker_id() { return omp_get_thread_num(); }
		int torc_i_num_workers() { return omp_get_max_threads(); }
		int torc_worker_id() { return omp_get_thread_num(); }
	#else
		int torc_i_worker_id() { return 0; }
		int torc_i_num_workers() { return 1; }
		int torc_worker_id() { return 0; }
	#endif

	#include <sys/time.h>
	double torc_gettime(){
    	struct timeval t;
    	gettimeofday(&t, NULL);
    	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
	}

#endif










/**********************************************/
/* Function call counters */
/**********************************************/
static pthread_mutex_t feval_m = PTHREAD_MUTEX_INITIALIZER;

static int l_nfeval = 0;
static int g_nfeval = 0;
static int t_nfeval = 0;

void inc_nfc()
{
    pthread_mutex_lock(&feval_m);
    l_nfeval++;
    pthread_mutex_unlock(&feval_m);
}

void reset_nfc_task()
{
    l_nfeval = 0;
}

void reset_nfc()
{
#if defined(_USE_TORC_)
    int i;

    for (i = 0; i < torc_num_nodes(); i++) {
        torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())reset_nfc_task, 0);
    }
    torc_waitall();
#else
    reset_nfc_task();
#endif
}

void get_nfc_task(int *x)
{
    *x = l_nfeval;
}

int get_nfc()
{
    int i;
    int c[1024]; /* MAX_NODES*/

#if defined(_USE_TORC_)
    for (i = 0; i < torc_num_nodes(); i++) {
        torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())get_nfc_task, 1,
                1, MPI_INT, CALL_BY_RES, &c[i]);
    }
    torc_waitall();
#else
    get_nfc_task(&c[0]);
#endif

    unsigned int s = 0;
    //printf("get_nfc:");
    for (i = 0; i < torc_num_nodes(); i++) {
		s += c[i];
    	// printf("+%d", c[i]);
    }
    g_nfeval = s;
    //printf("=%d\n", s);

    t_nfeval += g_nfeval;
    return g_nfeval;
}

int get_tfc()
{
    return t_nfeval;
}

/**********************************************/
/* Helper routines */
/**********************************************/
void print_matrix(char *title, double *v, int n)
{
    int i;

    /*    if (!display) return;*/

    printf("\n%s =\n\n", title);
    for (i = 0; i < n; i++) {
        printf("   %20.15lf\n", v[i]);
    }
    printf("\n");
}

void print_matrix_i(char *title, int *v, int n)
{
    int i;

    /*    if (!display) return;*/

    printf("\n%s =\n\n", title);
    for (i = 0; i < n; i++) {
        printf("  %8d\n", v[i]);
    }
    printf("\n");
}

void print_matrix_2d(char *title, double **v, int n1, int n2)
{
    int i, j;

    /*    if (!display) return;*/

    printf("\n%s =\n\n", title);
    for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
            printf("   %20.15lf", v[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void fprint_matrix_1d(FILE *fp, char *title, double *v, int n)
{
    int i;

    if (fp == stdout)
        fprintf(fp, "\n%s =\n\n", title);
    for (i = 0; i < n; i++) {
        fprintf(fp, "%12.4lf ", v[i]);
    }
    fprintf(fp, "\n");
}

void fprint_matrix_2d(FILE *fp, char *title, double **v, int n1, int n2)
{
    int i, j;

    if (fp == stdout)
        fprintf(fp, "\n%s =\n\n", title);
    for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
            fprintf(fp, "   %20.15lf", v[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

double compute_max(double *v, int n)
{
    int i;
    double vmax = v[0];
    for (i = 1; i < n; i++)
        if (v[i] > vmax) vmax = v[i];

    return vmax;
}

double compute_min(double *v, int n)
{
    int i;
    double vmin = v[0];
    for (i = 1; i < n; i++)
        if (v[i] < vmin) vmin = v[i];

    return vmin;
}

int compute_min_idx_i(int *v, int n)
{
    int i;
    double vmin = v[0];
    int idx = 0;

    for (i = 1; i < n; i++)
        if (v[i] < vmin) {
            vmin = v[i];
            idx = i;
        }

    return idx;
}

double compute_sum(double *v, int n)
{
    int i;
    double s = 0;
    for (i = 0; i < n; i++) s += v[i];

    return s;
}

double compute_mean(double *v, int n)
{
    int i;
    double s = 0;
    for (i = 0; i < n; i++) s += v[i];

    return s/n;
}

double compute_std(double *v, int n, double mean)
{
    int i;
    double s = 0;
    for (i = 0; i < n; i++) s += pow(v[i]-mean,2);

    return sqrt(s/(n-1));
}

/**********************************************/
/* Random number generators */
/**********************************************/
const gsl_rng_type *T;
gsl_rng **r;
int *local_seed;

void gsl_rand_init(int seed)
{
    int i, local_workers = torc_i_num_workers();
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = (gsl_rng **)malloc(local_workers*sizeof(gsl_rng *));
    for (i = 0; i < local_workers; i++) {
        r[i] = gsl_rng_alloc (T);
    }

    if (seed == 0) seed = time(0);

    for (i = 0; i < local_workers; i++) {
#if VERBOSE
        printf("node %d: initializing rng %d with seed %d\n", torc_node_id(), i, seed+i+local_workers*torc_node_id());
#endif
        gsl_rng_set(r[i], seed+i+local_workers*torc_node_id());
    }

    local_seed = (int *)malloc(local_workers*sizeof(int));
    for (i = 0; i < local_workers; i++) {
        local_seed[i] = seed+i+local_workers*torc_node_id();
    }
}

/* normal distribution random number N(mu,rho^2)*/
double normalrand(double mu, double var)
{
    double res;

    int me = torc_i_worker_id();
    res = mu + gsl_ran_gaussian(r[me], var);

    return res;

    /*    return mu + gsl_ran_gaussian(r, var);*/
}

/* uniform (flat) distribution random number between a and b */
double uniformrand(double a, double b)
{
    double res;

    int me = torc_i_worker_id();
    res = gsl_ran_flat(r[me], a, b);

    return res;

    /*    return gsl_ran_flat(r, a, b);*/
}

void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[])
{
    int me = torc_i_worker_id();
    gsl_ran_multinomial (r[me], K, N, q, nn);

    return;
}

void shuffle(int *perm, int N)
{
    int i;
    gsl_permutation * p = gsl_permutation_alloc (N);

    gsl_permutation_init (p);
#if VERBOSE
    gsl_permutation_fprintf (stdout, p, " %u");
    printf ("\n");
#endif
    int me = torc_i_worker_id();
    gsl_ran_shuffle (r[me], p->data, N, sizeof(size_t));
#if VERBOSE
    printf (" random permutation:");
    gsl_permutation_fprintf (stdout, p, " %u");
#endif
    for (i = 0; i <    N; i++)    perm[i] = p->data[i];
    gsl_permutation_free (p);
}

/**********************************************/
/* Multivariate Normal density function */
/**********************************************/
/*
 *  Multivariate Normal density function and random number generator
 *  Multivariate Student t density function and random number generator
 *  Wishart random number generator
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Copyright (C) 2006  Ralph dos Santos Silva
 */
int mvnrnd_silva(const gsl_rng *r, const gsl_vector *mean, const gsl_matrix *var, gsl_vector *result)
{
    size_t n = mean->size;
    /* multivariate normal distribution random number generator */
    /*
     *    n    dimension of the random vetor
     *    mean    vector of means of size n
     *    var     variance matrix of dimension n x n
     *    result  output variable with a sigle random vector normal distribution generation
     */

    unsigned int k;
    gsl_matrix *work = gsl_matrix_alloc(n,n);

    gsl_matrix_memcpy(work,var);
    /*    printf("I am here\n");*/
    gsl_linalg_cholesky_decomp(work);

    for(k=0; k<n; k++)
        gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result);
    gsl_vector_add(result,mean);

    gsl_matrix_free(work);

    return 0;
}

/*
 *  @title multivariate normal random variables
 *  @author Carl Boettiger, <cboettig@gmail.com>
 *
 *  Based on the R function rmvnorm, from the mvtnorm package
 *  by Friedrich Leisch and Fabian Scheipl, implemented
 *  using the GSL libraries
 */

int mvnrnd_cboet(gsl_rng * rng, const gsl_vector * mean, gsl_matrix * covar, gsl_vector * ANS)
{
    unsigned int i;
    size_t n = mean->size;

    /* Calculate eigenvalues and eigenvectors of covar matrix */
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (n);
    gsl_eigen_symmv (covar, eval, evec, w);
    gsl_eigen_symmv_free (w);
    /*    gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);*/



    /* Setup for: evec * matrix(diag(eval)) * transpose(evec)  */
    gsl_matrix *eval_mx = gsl_matrix_calloc (n, n);
    gsl_matrix * x_M = gsl_matrix_alloc (n,n);
    gsl_matrix * x_M_x = gsl_matrix_alloc (n,n);


    gsl_vector_view diagonal = gsl_matrix_diagonal(eval_mx);
    gsl_vector_memcpy(&diagonal.vector, eval);
    for(i=0;i<n;i++)
    {
        gsl_vector_set( &diagonal.vector,
                i,
                sqrt( gsl_vector_get(&diagonal.vector, i) )
                );
    }



    /* evec * matrix(diag(eval)) * transpose(evec)  */
    /*    gsl_blas_dsymm (CblasLeft, CblasUpper, 1.0, evec, eval_mx, 0.0, x_M);*/

    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
            1.0, evec, eval_mx, 0.0, x_M);
    gsl_blas_dgemm (CblasNoTrans, CblasTrans,
            1.0, x_M, evec, 0.0, x_M_x);


    gsl_matrix_free(x_M);
    gsl_matrix_free(eval_mx);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);

    gsl_vector * rnorms = gsl_vector_alloc(n);
    for(i=0;i<n;i++)
    {
        gsl_vector_set
            ( rnorms, i,
              gsl_ran_gaussian_ziggurat(rng, 1)
            );
    }

    gsl_blas_dgemv( CblasTrans, 1.0, x_M_x, rnorms, 0, ANS);
    gsl_vector_add(ANS, mean);
    gsl_matrix_free(x_M_x);
    gsl_vector_free(rnorms);

    return 0;
    /* answer provided through pass by reference */
}

int mvnrnd_gsl(gsl_rng * rng, const gsl_vector * mean, gsl_matrix * covar, gsl_vector * ANS)
{
#if 1
    return mvnrnd_silva(rng, mean, covar, ANS);
#else
    return mvnrnd_cboet(rng, mean, covar, ANS);
#endif
}

int mvnrnd(double *mean, double *sigma, double *out, int N)
{
    int res;

    gsl_vector_view mean_view = gsl_vector_view_array(mean, N);
    gsl_matrix_view sigma_view = gsl_matrix_view_array(sigma, N,N);
    gsl_vector_view out_view = gsl_vector_view_array(out, N);

    int me = torc_i_worker_id();
    res = mvnrnd_gsl(r[me], &mean_view.vector, &sigma_view.matrix, &out_view.vector);

    return res;
}

#if 1
void aux_init()
{
#if defined(_USE_TORC_)
    torc_register_task((void *)reset_nfc_task);
    torc_register_task((void *)get_nfc_task);
#endif
}
#endif

/* Copyright (C) 2006  Ralph dos Santos Silva */
/* http://lists.gnu.org/archive/html/help-gsl/2006-04/txtdb8Hdlx9uA.txt */

static double dmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var, int lognorm)
{
    /* multivariate normal density function    */
    /*
     *    n    dimension of the random vetor
     *    mean    vector of means of size n
     *    var    variance matrix of dimension n x n
     */

    int s;
    double ax,ay;
    gsl_vector *ym, *xm;
    gsl_matrix *work = gsl_matrix_alloc(n,n), *winv = gsl_matrix_alloc(n,n);
    gsl_permutation *p = gsl_permutation_alloc(n);

    gsl_matrix_memcpy( work, var );
    gsl_linalg_LU_decomp( work, p, &s );
    gsl_linalg_LU_invert( work, p, winv );
    ax = gsl_linalg_LU_det( work, s );
    gsl_matrix_free( work );
    gsl_permutation_free( p );

    xm = gsl_vector_alloc(n);
    gsl_vector_memcpy( xm, x);
    gsl_vector_sub( xm, mean );
    ym = gsl_vector_alloc(n);
    gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
    gsl_matrix_free( winv );
    gsl_blas_ddot( xm, ym, &ay);
    gsl_vector_free(xm);
    gsl_vector_free(ym);

    if (!lognorm)
        ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );
    else
        ay = 0.5*(-ay - n*log(2*M_PI) - log(ax));    /* PA */

    return ay;
}

static double gsl_dmvnorm(int n, double *xv, double *mv, double *vm, int lognorm)
{
    double result;
    int i, j;

    gsl_vector *x = gsl_vector_calloc(n), *mean = gsl_vector_calloc(n);
    gsl_matrix *var = gsl_matrix_calloc(n,n);

    for (i = 0; i < n; i++)
        gsl_vector_set(x,i,xv[i]);

    if (mv != NULL) {
        for (i = 0; i < n; i++)
            gsl_vector_set(mean,i,mv[i]);
    }

    if (vm != NULL) {
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                gsl_matrix_set(var,i,j,vm[i*n+j]);
    }
    else {
        for (i = 0; i < n; i++)
            gsl_matrix_set(var,i,i, 1.0);
    }

#if 0
    gsl_vector_fprintf(stdout, x, "%f");
    gsl_vector_fprintf(stdout, mean, "%f");
    gsl_matrix_fprintf(stdout, var, "%f");
#endif

    result = dmvnorm(n,x,mean,var,lognorm);

    gsl_vector_free(x);
    gsl_vector_free(mean);
    gsl_matrix_free(var);

    return result;

}

double mvnpdf(int n, double *xv, double *mv, double *vm)
{
    double result = gsl_dmvnorm(n, xv, mv, vm, 0);
    return result;
}

double logmvnpdf(int n, double *xv, double *mv, double *vm)
{
    double result = gsl_dmvnorm(n, xv, mv, vm, 1);
    return result;
}


#include "thirdparty/truncated_normal.c"

/******************************************************************************/
double truncated_normal_pdf (double x, double mu, double sigma, double a, double b)
{
    double pdf;

    pdf = truncated_normal_ab_pdf(x, mu, sigma, a, b);

    return pdf;
}

/******************************************************************************/
double truncated_normal_rand (double mu, double sigma, double a, double b)
{
    int me = torc_i_worker_id();

    double x = truncated_normal_ab_sample(mu, sigma, a, b, &local_seed[me]);

    return x;
}

/******************************************************************************/
double truncated_lognormal_pdf (double x, double mu, double sigma, double a, double b)
{
    double pdf;

    pdf = log_normal_truncated_ab_pdf(x, mu, sigma, a, b);

    return pdf;
}
















//----------------------------------------------------------------------------------
//
//
void call_gsl_rand_init()
{
    // printf("CALLING gsl_rand_init() on node %d\n", torc_node_id()); fflush(0);
    gsl_rand_init(data.seed);
}



void spmd_gsl_rand_init()
{
	#if defined(_USE_TORC_)
    	int i;
    	for (i = 0; i < torc_num_nodes(); i++) {
        	torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())call_gsl_rand_init, 0);
    	}
    	torc_waitall();
	#else
		call_gsl_rand_init();
	#endif
}



void call_print_matrix_2d()
{
    printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
    print_matrix_2d((char *)"runinfo.SS", runinfo.SS, data.Nth, data.Nth);
    printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
}



void spmd_print_matrix_2d()
{
	#if defined(_USE_TORC_)
		int i;
		for (i = 0; i < torc_num_nodes(); i++) {
			torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())call_print_matrix_2d, 0);
		}
		torc_waitall();
	#else
		call_print_matrix_2d();
	#endif
}



void call_update_gdata()    /* step for p[j]*/
{
	#if defined(_USE_TORC_)
		MPI_Bcast(runinfo.SS[0], data.Nth*data.Nth, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(runinfo.p, data.MaxStages, MPI_DOUBLE, 0, MPI_COMM_WORLD);    /* just p[Gen]*/
		MPI_Bcast(&runinfo.Gen, 1, MPI_INT, 0, MPI_COMM_WORLD);
	#endif
}



void spmd_update_gdata()    /* step*/
{
	#if defined(_USE_TORC_)
		int i;
		if (torc_num_nodes() == 1) return;
		for (i = 0; i < torc_num_nodes(); i++) {
			torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())call_update_gdata, 0);
		}
		torc_waitall();
	#endif
}









