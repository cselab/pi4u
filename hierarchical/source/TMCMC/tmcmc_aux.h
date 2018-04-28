/*
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#ifndef _TMCMC_AUX_H_
#define _TMCMC_AUX_H_



void multinomialrand(size_t K, unsigned int N, double q[], unsigned int nn[]);
int mvnrnd(double *mean, double *sigma, double *out, int N);


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


/*** AUX ***/
void inc_nfc();
void get_nfc_task(int *);
int get_nfc();
void reset_nfc_task();
void reset_nfc();
int get_tfc();

void call_gsl_rand_init();
void spmd_gsl_rand_init();
void call_print_matrix_2d();
void spmd_print_matrix_2d();
void call_update_gdata();
void spmd_update_gdata();




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

	int torc_node_id();
	int torc_num_nodes();

	#if defined(_USE_OPENMP_)
		int torc_i_worker_id();
		int torc_i_num_workers();
		int torc_worker_id();
	#else
		int torc_i_worker_id();
		int torc_i_num_workers();
		int torc_worker_id();
	#endif

	double torc_gettime();

#endif


#include <gsl/gsl_rng.h>
extern gsl_rng   				**r;
extern int   					*local_seed;




#endif
