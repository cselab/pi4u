/*
 *  tmcmc_kriging.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <time.h>

#include "engine_tmcmc.h"
#include "gsl_headers.h"
#include "spawner.h"


#ifndef MAX
#define MAX(a,b) ((a)>(b)) ? (a):(b)
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)) ? (a):(b)
#endif

// this must be global information (MPI)
static const double ERROR_THRESHOLD=0.1;	// set once
static const int PREVIOUS_STAGES=2;

static double v_min=0, v_max=0;

/*static*/ void call_update_nn_gdata()
{
	MPI_Bcast(&v_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&v_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

static void spmd_update_nn_gdata()
{
	if (torc_num_nodes() == 1) return;

	int i;
	for (i = 0; i < torc_num_nodes(); i++)
		torc_create_ex(i*torc_i_num_workers(), 1, call_update_nn_gdata, 0);
	torc_waitall();
}


int surrogate_is_good_estimate(double estimate, double error)
{
	// the avg_error might be high only because of some extreme values
	// should we look at the median + some classification of estimated values (in addition to quantiles) ?

	if ((error < ERROR_THRESHOLD) && (v_min < estimate) && (estimate < v_max))
		return 1;

	return 0;
}

void compute_surrogate_bounds()
{
	int N = full_db.entries;
	int F = 1;

	double *T = calloc(1, N*F*sizeof(double));

	int i;
	
	int real_i = 0, real_N;
	for (i = 0; i < N; i++)
	{
		if (i < N-PREVIOUS_STAGES*data.PopSize) continue;	// how many past stages should we use? 
		if (full_db.entry[i].surrogate == 1) continue;

		T[real_i] = full_db.entry[i].F;
		real_i++;
	}
	real_N = real_i;

	// compute percentile - use only real model evaluations
	double * loglik_array = T;
	int count = real_N; 

	gsl_sort(loglik_array, 1, count);
	double p_001 = gsl_stats_quantile_from_sorted_data(loglik_array, 1, count, 0.01);
	double p_005 = gsl_stats_quantile_from_sorted_data(loglik_array, 1, count, 0.05);
	double p_095 = gsl_stats_quantile_from_sorted_data(loglik_array, 1, count, 0.95);
	double p_099 = gsl_stats_quantile_from_sorted_data(loglik_array, 1, count, 0.99);
	double loglik_min = loglik_array[0];
	double loglik_max = loglik_array[count-1];

//	v_min = loglik_min - 0.01*(loglik_max-loglik_min);
//	v_max = loglik_max + 0.01*(loglik_max-loglik_min);
	v_min = p_001; // if the surrogate's value is too big/small, perhaps it's an oscillation
	v_max = p_099;

	free(T);

	spmd_update_nn_gdata();

	printf("STAT: v_min = %f v_max = %f\n", v_min, v_max);
}

// Build a surrogate each time the function is called (dumping files is more expensive)
void get_surrogate_estimate(int gen_id, int chain_id, int task_id, double *x,
	double *pred_y, double *err, double *leader)
{
	int i;
	double t;

	int D = data.Nth;

//	t = clock();

	char command[1024];
	sprintf(command, "python scripts/kriging_doall.py \'%d %d %d %d ", gen_id, chain_id, task_id, D);
	for (i=0; i<D; i++)
		sprintf(command, "%s %lf", command, leader[i]);
	for (i=0; i<D; i++)
		sprintf(command, "%s %lf", command, x[i]);
	sprintf(command, "%s\' 2>/dev/null\n", command);

	FILE * pipe = popen(command, "r");

	if (!pipe)
	{
	    printf("Cannot popen %s\n", command);
		return;
	}

	t = -1;
	*pred_y = -1;
	*err = 1e12;

	char buffer[512];
	while(fgets(buffer, 512, pipe))
	{
//		printf("(in C from python) %s\n", buffer);
		char word[512];
		sscanf(buffer, "%s ", word);
		if (strcmp(word, "Prediction:") == 0)
			sscanf(buffer, "%*s %lf", pred_y);
		if (strcmp(word, "Error:") == 0)
			sscanf(buffer, "%*s %lf", err);
		if (strcmp(word, "Elapsed") == 0)
			sscanf(buffer, "%*s %*s %lf", &t);
	}

	pclose(pipe);
		
//	t = (clock() - t)/CLOCKS_PER_SEC;
	printf("Prediction elapsed time: %.2lf sec\n", t);
}

void destroy_net()
{
}


