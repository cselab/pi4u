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

// the surrogate is the same for the whole chain
void surrogate_build(int gen_id, int chain_id, double *x)
{
	int D = data.Nth;

	char filename1[512], filename2[512];
	FILE *fp;
	int i;

	sprintf(filename1, "leader_%03d_%03d.txt", gen_id, chain_id);
	fp = fopen(filename1, "w");
	for (i=0; i<D; i++)
		fprintf(fp, "%lf ", x[i]);
	fprintf(fp, "\n");
	fclose(fp);

	// Call a script to train the kriging model
	char line[1024];
	char *largv[64];

	int rf = fork();
	if (rf < 0)
	{
		printf("spawner(%d, %d): fork failed!!!!\n", gen_id, chain_id);
		fflush(0);
	}
	if (rf == 0)
	{
		strcpy(line, "");
		sprintf(line, "python scripts/kriging_train.py %d %d\n", gen_id, chain_id);

		sprintf(filename2, "train_output_%03d_%03d.txt", gen_id, chain_id);

		int fd = open(filename2, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient

		parse(line, largv);	/* prepare argv */
		execvp(*largv, largv);
	}
	waitpid(rf, NULL, 0);
	sync();

	rmrf(filename1);
	rmrf(filename2);
}

void get_surrogate_estimate(int gen_id, int chain_id, int task_id, double *x, int n,
	double *pred_y, double *err)
{
	char filename1[512], filename2[512], filename3[512];
	FILE * fp;
	int i;
	
	// write x to file
	sprintf(filename1, "kriging_pred_x_%03d_%03d_%03d.txt", gen_id, chain_id, task_id);
	fp = fopen(filename1, "w");
	for(i=0; i<n; i++)
		fprintf(fp,"%lf ", x[i]);
	fprintf(fp,"\n");
	fclose(fp);

	// Obtain the prediction and the error
	char line[1024];
	char *largv[64];

	int rf = fork();
	if (rf < 0)
	{
		printf("spawner(%d, %d): fork failed!!!!\n", gen_id, chain_id);
		fflush(0);
	}
	if (rf == 0)
	{
		strcpy(line, "");
		sprintf(line, "python scripts/kriging_pred.py %d %d %d\n", gen_id, chain_id, task_id);

		sprintf(filename2, "pred_output_%03d_%03d_%03d.txt", gen_id, chain_id, task_id);

		int fd = open(filename2, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient

		parse(line, largv);	/* prepare argv */
		execvp(*largv, largv);
	}
	waitpid(rf, NULL, 0);
	sync();

	// Read the prediction and the error from the file
	double krig_pred_y, krig_err;
	sprintf(filename3, "kriging_pred_y_%03d_%03d_%03d.txt", gen_id, chain_id, task_id);
	fp = fopen(filename3, "r");
	if (!fp)
	{
		printf("Error in prediction for %s\n", filename3);
	
		*pred_y = -1;
		*err = 1e12;
	}
	else
	{
		fscanf(fp, "%lf\n%lf\n", &krig_pred_y, &krig_err);
		fclose(fp);
	
		printf("evaluate_surrogate -> %f with error %f\n", krig_pred_y, krig_err);
	
		*pred_y = krig_pred_y;
		*err = krig_err;

		rmrf(filename3);
	}

	rmrf(filename1);
	rmrf(filename2);
}

void destroy_net()
{
}


