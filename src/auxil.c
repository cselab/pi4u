/*
 *  auxil.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
//#include <torc.h>

#include <math.h>
#include <time.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>


/**********************************************/
// Function call counters
/**********************************************/

#ifdef USE_TORC
static pthread_mutex_t feval_m = PTHREAD_MUTEX_INITIALIZER;
#endif

static int l_nfeval = 0;
static int g_nfeval = 0;
static int t_nfeval = 0;

void inc_nfc()
{
#ifdef USE_TORC
	pthread_mutex_lock(&feval_m);
#endif
	l_nfeval++;
#ifdef USE_TORC
	pthread_mutex_unlock(&feval_m);
#endif
}

void reset_nfc_task()
{
	l_nfeval = 0;
}

void reset_nfc()
{
#ifdef USE_TORC
	int i;

	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, (void *)reset_nfc_task, 0);
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
#ifdef USE_TORC
	int i;
	int c[64]; // MAX_NODES

	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, (void *)get_nfc_task, 1,
		1, MPI_INT, CALL_BY_RES, &c[i]);
	}
	torc_waitall();

	unsigned int s = 0;
	printf("get_nfc:");
	for (i = 0; i < torc_num_nodes(); i++) {
		s += c[i];
		printf("+%d", c[i]);
	}
	g_nfeval = s;
	printf("=%d\n", s);
#else
	g_nfeval = l_nfeval;
#endif
	t_nfeval += g_nfeval;
	return g_nfeval;
}

int get_tfc()
{
	return t_nfeval;
}

/**********************************************/
// Helper routines
/**********************************************/

void print_matrix(char *title, double *v, int n)
{
	int i;

//	if (!display) return;

	printf("\n%s =\n\n", title);
	for (i = 0; i < n; i++) {
		printf("   %20.15lf\n", v[i]);
	}
	printf("\n");
}

void print_matrix_i(char *title, int *v, int n)
{
	int i;

//	if (!display) return;

	printf("\n%s =\n\n", title);
	for (i = 0; i < n; i++) {
		printf("  %8d\n", v[i]);
	}
	printf("\n");
}

void print_matrix_2d(char *title, double **v, int n1, int n2)
{
	int i, j;

//	if (!display) return;

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
// Random number generators
/**********************************************/

const gsl_rng_type *T;
gsl_rng *r;
#ifdef USE_TORC
static pthread_mutex_t _rm = PTHREAD_MUTEX_INITIALIZER;
#endif
static int rng_initialized = 0;

void gsl_rand_init(int seed)
{
#ifdef USE_TORC
		pthread_mutex_lock(&_rm);
#endif
		if (rng_initialized == 0) {
			rng_initialized = 1;
		}
		else {
#ifdef USE_TORC
			pthread_mutex_unlock(&_rm);
#endif
			return;
		}

		gsl_rng_env_setup();

		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
		if (seed == 0) 
#ifdef USE_TORC
			gsl_rng_set(r, time(0)+10*torc_node_id());
#else
			gsl_rng_set(r, time(0));
#endif
		else
			gsl_rng_set(r, seed);

#ifdef USE_TORC
		pthread_mutex_unlock(&_rm);
#endif
}

/* normal distribution random number N(mu,rho^2)*/
double normalrand(double mu, double var)
{
	double res;

#ifdef USE_TORC
	pthread_mutex_lock(&_rm);
#endif
	res = mu + gsl_ran_gaussian(r, var);
#ifdef USE_TORC
	pthread_mutex_unlock(&_rm);
#endif
	return res;

//	return mu + gsl_ran_gaussian(r, var);
}

/* uniform (flat) distribution random number between a and b */
double uniformrand(double a, double b)
{
	double res;

#ifdef USE_TORC
	pthread_mutex_lock(&_rm);
#endif
	res = gsl_ran_flat(r, a, b);
#ifdef USE_TORC
	pthread_mutex_unlock(&_rm);
#endif
	return res;

//	return gsl_ran_flat(r, a, b);
}


void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[])
{
#ifdef USE_TORC
	pthread_mutex_lock(&_rm);
#endif
	gsl_ran_multinomial (r, K, N, q, nn);
#ifdef USE_TORC
	pthread_mutex_unlock(&_rm);
#endif
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
	gsl_ran_shuffle (r, p->data, N, sizeof(size_t));
#if VERBOSE
	printf (" random permutation:");
	gsl_permutation_fprintf (stdout, p, " %u");
#endif
	for (i = 0; i <	N; i++)	perm[i] = p->data[i];
	gsl_permutation_free (p);
}


/**********************************************/
// Multivariate Normal density function
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
	*	n	dimension of the random vetor
	*	mean    vector of means of size n
	*	var     variance matrix of dimension n x n
	*	result  output variable with a sigle random vector normal distribution generation
	*/

	unsigned int k;
	gsl_matrix *work = gsl_matrix_alloc(n,n);

	gsl_matrix_memcpy(work,var);
//	printf("I am here\n");
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
//	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);



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
//	gsl_blas_dsymm (CblasLeft, CblasUpper, 
//					1.0, evec, eval_mx, 0.0, x_M);

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

#ifdef USE_TORC
	pthread_mutex_lock(&_rm);
#endif
	res = mvnrnd_gsl(r, &mean_view.vector, &sigma_view.matrix, &out_view.vector);
#ifdef USE_TORC
	pthread_mutex_unlock(&_rm);
#endif

	return res;
}


#if 1
void aux_init()
{
#ifdef USE_TORC
	torc_register_task(reset_nfc_task);
	torc_register_task(get_nfc_task);
#endif
}
#endif
