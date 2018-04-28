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

#include "../priors/myrand.h"







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



void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[])
{
    int me = torc_i_worker_id();
    gsl_ran_multinomial (r[me], K, N, q, nn);

    return;
}






int mvnrnd(double *mean, double *sigma, double *out, int N)
{
    int res;

    gsl_vector_view mean_view 	= gsl_vector_view_array(mean, N);
    gsl_matrix_view sigma_view 	= gsl_matrix_view_array(sigma, N,N);
    gsl_vector_view out_view 	= gsl_vector_view_array(out, N);

    int me = torc_i_worker_id();

    gsl_matrix *L = gsl_matrix_alloc(N,N);
    gsl_matrix_memcpy( L, &sigma_view.matrix);
    gsl_linalg_cholesky_decomp( L );


	res = gsl_ran_multivariate_gaussian( r[me], &mean_view.vector, L, &out_view.vector);

    return res;
}







// TODO: move this to priors
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
//		RNG initialization
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







//----------------------------------------------------------------------------------
//				Print functions
//		
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





