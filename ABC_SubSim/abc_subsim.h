/*
 *  abc_subsim.h
 *  Pi4U
 *
 *  Created by Lina Kulakova on 1/1/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>

class abc_subsim
{
private:
	int population_size;
	std::string posterior_filename;
	std::string restart_filename;
	double tolerance_final;
	double acc_rate;
	int max_level;
	double weights_std;
	double tolerance_decrease;
	int seed;
	int pivot;
	int chain_length;
	int chain_length_test;
	int restart;

public:
	abc_subsim(int population_size, const char * posterior_filename, double tolerance_final, int max_level, double acc_rate,
			double weights_std,	double tolerance_decrease, int seed, int restart, const char * restart_filename):
		population_size(population_size), posterior_filename(posterior_filename), tolerance_final(tolerance_final),
		max_level(max_level), acc_rate(acc_rate), weights_std(weights_std), tolerance_decrease(tolerance_decrease),
		seed(seed), restart(restart), restart_filename(restart_filename)
	{
		chain_length = 5; // keep 1/chain_length best individuals		// user input
		chain_length_test = 5; //std::min(std::max(chain_length,10),10);		// user input
		pivot = floor(population_size/chain_length);
	}

	void run();

	// peh: function address should be visible across different processes
	static
	void run_mh(int *p_chain_length,
			const double * seed,
			double *p_discrepancy_seed,
			double *p_tolerance,
			double * acc_rate_array_val,
			int *p_seed_ind, int *p_level);

	static
	void run_mh_test(int chain_length,
			const double * seed_v,
			double discrepancy_seed,

			const gsl_matrix * covariance,
			double tolerance,

			double * acc_rate_array_val,
			int seed_ind,
			int level);

	void read_restart_file(std::vector< gsl_vector * > & sample_array, gsl_vector * discrepancy_array);

	void sample_from_prior_in_parallel(std::vector< gsl_vector * > & sample_array, gsl_vector * discrepancy_array);

	void tune_covariance(std::vector< gsl_vector * > const& sample_array_prev,
			const gsl_vector * discrepancy_array_prev, double tolerance_curr, std::vector< gsl_vector * > & sample_array_curr,
			gsl_vector * discrepancy_array_curr, gsl_vector * acc_rate_array,
			gsl_matrix * covariance, const size_t * perm, int level);

	static
	void tune_covariance_task(
			int *pchain_length_test,
			const double *sample_v,
			double *pdiscrepancy_v,
			double *ptolerance_curr,
			double *acc_rate_array_val,
			double *pscale,
			int *pnumber,
			int *plevel);

	void tune_covariance_in_parallel(
			std::vector< gsl_vector * > const& sample_array_prev,
			const gsl_vector * discrepancy_array_prev,
			double tolerance_curr,
			std::vector< gsl_vector * > & sample_array_curr,
			gsl_vector * discrepancy_array_curr,
			const size_t * perm,
			int level);

	void tune_covariance_in_parallel_v0(
			std::vector< gsl_vector * > const& sample_array_prev,
			const gsl_vector * discrepancy_array_prev,
			double tolerance_curr,
			std::vector< gsl_vector * > & sample_array_curr,
			gsl_vector * discrepancy_array_curr,
			const size_t * perm,
			int level);

	void run_mh_in_parallel(int chain_length,
				const size_t * perm,
				std::vector< gsl_vector * > & sample_array_prev,
				gsl_vector * discrepancy_array_prev,
				double tolerance_curr,
				gsl_vector * acc_rate_array,
				int level);

	void write_array_to_file(std::vector< gsl_vector * > const& sample, const gsl_vector * discrepancy, int level);
};

