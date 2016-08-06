/*
 *  random.h
 *  Pi4U
 *
 *  Created by Lina Kulakova on 1/1/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#if defined(_USE_TORC_)
#include <mpi.h>
extern "C"
{
#include <torc.h>
}
#else
int torc_i_worker_id()
{
	return 0;
}
#endif


/*
 * peh: (I) class myrandom
 */

#define MAX_WORKERS	48

class myrandom {
public:
	static myrandom *me;
	gsl_rng ** r;

	myrandom() {
		r = (gsl_rng **)malloc(MAX_WORKERS*sizeof(gsl_rng *));
		for (int i = 0; i < MAX_WORKERS; i++)
			r[i] = gsl_rng_alloc(gsl_rng_ranlux389);
		myrandom::me = this;
	}

	void set(int seed)
	{
		for (int i = 0; i < MAX_WORKERS; i++)
			gsl_rng_set(r[i], seed+i);
	}

	static void set_task(int *pseed)
	{
		int seed = *pseed;
		if (seed == 0) seed = time(0);
		//printf("CALLING set_task(%d) on node %d\n", seed, torc_node_id()); fflush(0);
		for (int i = 0; i < MAX_WORKERS; i++)
			myrandom::me->set(seed+MAX_WORKERS*torc_node_id()+i);
	}

	void spmd_set(int seed)
	{
#if defined(_USE_TORC_)
		for (int i = 0; i < torc_num_nodes(); i++) {
			torc_create_ex(i*torc_i_num_workers(), 1, (void (*)())myrandom::set_task, 1,
					1, MPI_INT, CALL_BY_COP,
					&seed);
		}
		torc_waitall();
#else
		myrandom::set_task(&seed);
#endif
	}
};

myrandom * myrandom::me = NULL;	// tricky but not actually necessary as the rng object is global

myrandom rng;

double gsl_ran_flat_logpdf(double x, double a, double b, int * isinf)
{
	double res = 0.0;

	if(x < a || b < x)
	{
		res = -1;
		*isinf = 1;
	}
	else
	{
		res = log(gsl_ran_flat_pdf(x,a,b));
		*isinf = 0;
	}

	return res;
}

// code from http://lists.gnu.org/archive/html/help-gsl/2006-04/txtdb8Hdlx9uA.txt
int ran_mvn(const gsl_vector * mean, const gsl_matrix * covariance, gsl_vector * sample)
{
	int n = mean->size;
	gsl_matrix * chol = gsl_matrix_alloc(n,n);

	gsl_matrix_memcpy(chol, covariance);
	gsl_linalg_cholesky_decomp(chol);

	// GSL is thread-safe
	for(int k=0; k<n; k++)
		gsl_vector_set(sample, k, gsl_ran_ugaussian(rng.r[torc_i_worker_id()]));

	gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, chol, sample);
	gsl_vector_add(sample, mean);

	gsl_matrix_free(chol);

	return 0;
}

// code from http://lists.gnu.org/archive/html/help-gsl/2006-04/txtdb8Hdlx9uA.txt
double ran_mvn_pdf(const gsl_vector * x, const gsl_vector * mean, const gsl_matrix * covariance)
{
	int n = mean->size;

	// workspace
	gsl_matrix * Sigma = gsl_matrix_alloc(n,n);
	gsl_matrix * Sigma_inv = gsl_matrix_alloc(n,n);
	gsl_permutation * p = gsl_permutation_alloc(n);
	gsl_vector * xm = gsl_vector_alloc(n);
	gsl_vector * ym = gsl_vector_alloc(n);

	int s;
	gsl_matrix_memcpy(Sigma, covariance);
	gsl_linalg_LU_decomp(Sigma, p, &s);
	gsl_linalg_LU_invert(Sigma, p, Sigma_inv);
	double det = gsl_linalg_LU_det(Sigma, s);

	gsl_vector_memcpy(xm, x);
	gsl_vector_sub(xm, mean);

	gsl_blas_dsymv(CblasUpper, 1.0, Sigma_inv, xm, 0.0, ym);

	double sos;
	gsl_blas_ddot(xm, ym, &sos);

	double res = exp(-0.5*sos)/sqrt(pow((2*M_PI),n)*det);

	// free workspace
	gsl_matrix_free(Sigma);
	gsl_matrix_free(Sigma_inv);
	gsl_permutation_free(p);
	gsl_vector_free(xm);
	gsl_vector_free(ym);

	return res;
}

double ran_mvn_logpdf(const gsl_vector * x, const gsl_vector * mean, const gsl_matrix * covariance)
{
	int n = mean->size;

	// workspace
	gsl_matrix * Sigma = gsl_matrix_alloc(n,n);
	gsl_matrix * Sigma_inv = gsl_matrix_alloc(n,n);
	gsl_permutation * p = gsl_permutation_alloc(n);
	gsl_vector * xm = gsl_vector_alloc(n);
	gsl_vector * ym = gsl_vector_alloc(n);

	int s;
	gsl_matrix_memcpy(Sigma, covariance);
	gsl_linalg_LU_decomp(Sigma, p, &s);
	gsl_linalg_LU_invert(Sigma, p, Sigma_inv);
	double det = gsl_linalg_LU_det(Sigma, s);

	gsl_vector_memcpy(xm, x);
	gsl_vector_sub(xm, mean);

	gsl_blas_dsymv(CblasUpper, 1.0, Sigma_inv, xm, 0.0, ym);

	double sos;
	gsl_blas_ddot(xm, ym, &sos);

	double res = -0.5*( sos + n*log(2*M_PI) + log(det) );

	// free workspace
	gsl_matrix_free(Sigma);
	gsl_matrix_free(Sigma_inv);
	gsl_permutation_free(p);
	gsl_vector_free(xm);
	gsl_vector_free(ym);

	return res;
}
