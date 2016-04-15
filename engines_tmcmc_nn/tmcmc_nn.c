/*
 *  tmcmc_nn.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include "engine_tmcmc.h"
#include "gsl_headers.h"

#if 0
#include "nnet.h"
#else
void build_network(char *netname, double *Xin, double *Tin, int D, int F, int N, double *Eout);	// #layers, #neurons, other parameters?
double predict_network(char *netname, double *Xin, int D);
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b)) ? (a):(b)
#endif

#ifndef MIN
#define MIN(a,b) ((a)<(b)) ? (a):(b)
#endif

// this must be global information (MPI)
static const double ERROR_THRESHOLD=5;	// set once

static double avg_error=100, p_005= 0, p_095=0;

/*static*/ void call_update_nn_gdata()
{
        MPI_Bcast(&avg_error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&p_005, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&p_095, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

static void spmd_update_nn_gdata()
{
        int i;
	if (torc_num_nodes() == 1) return;
        for (i = 0; i < torc_num_nodes(); i++) {
                torc_create_ex(i*torc_i_num_workers(), 1, call_update_nn_gdata, 0);
        }
	torc_waitall();
}


int surrogate_is_good_estimate(double estimate)
{
	// the avg_error might be high only because of some extreme values
	// should we look at the median + some classification of estimated values (in addition to quantiles) ?
	if ((avg_error < ERROR_THRESHOLD) && (p_005 < estimate) && (estimate < p_095))
		return 1;

	return 0;
}

void surrogate_build(int Gen)
{
	int D = data.Nth;
	int N = full_db.entries;
	int F = 1;
	int previous_stages = 2;

	double *X = calloc(1, N*D*sizeof(double));
	double *T = calloc(1, N*F*sizeof(double));
	double *E = calloc(1, N*F*sizeof(double));
	int *idx  = calloc(1, N*sizeof(int));

	int real_i = 0;
	int i;
	for (i = 0; i < N; i++)
	{
		double x[D], t;

		idx[i] = -1;

		int j;
		for (j = 0; j < D; j++) x[j] = full_db.entry[i].point[j];
		t = full_db.entry[i].F;

		if (i < N-previous_stages*data.PopSize) continue;	// how many past stages should we use? 
		if (full_db.entry[i].surrogate == 1) continue;

		for (j = 0; j < D; j++) X[real_i*D+j] = x[j];
		T[real_i] = t;
		idx[i] = real_i;

		real_i++;
	}
	int real_N = real_i;

	char netname[256];
	sprintf(netname, "myNN_%03d.txt", Gen);
	build_network(netname, X, T, D, F, real_N, E);
//	build_network("myNN.txt", X, T, D, F, real_N, E);

        double max_re = 0;
        double sum_re = 0;

	#define NBINS   20
	int bins[NBINS+1];
	for (i = 0; i < NBINS+1; i++) bins[i]=0;

	for (i = 0; i < N; i++) {
		if (idx[i] == -1) continue;

		double error = E[idx[i]];

		if (max_re > error) max_re = error;
		sum_re += error;

		int cbin = (int)error/(100.0/NBINS);
		if (cbin > NBINS) cbin = NBINS;
		if (error > 100) cbin = NBINS;
		bins[cbin]++;
        }

	double avg_re = sum_re/real_N;

	printf("max abs rel error = %.2f%%\n", max_re);
	printf("avg abs rel error = %.2f%%\n", avg_re);

	int sbin = 0;
	for (i = 0; i < NBINS+1; i++) sbin+=bins[i];
	for (i = 0; i < NBINS; i++)
	{
		printf("[%3.0f - %3.0f) -> %d (%.2f%%)\n", (100.0/NBINS)*i, (100.0/NBINS)*(i+1), bins[i], 100.0*bins[i]/sbin);
	}
	printf("[%3.0f - inf) -> %d (%.2f%%)\n", 100.0, bins[NBINS], 100.0*bins[NBINS]/sbin);

	{
	// compute percentile - use only real model evaluations
	double * loglik_array = T;
	int count = real_N; 

	gsl_sort(loglik_array, 1, count);
	p_005 = gsl_stats_quantile_from_sorted_data(loglik_array, 1, count, 0.05);
	p_095 = gsl_stats_quantile_from_sorted_data(loglik_array, 1, count, 0.95);
	}

	free(X);
	free(T);
	free(E);
	free(idx);


	avg_error = avg_re;

	spmd_update_nn_gdata();

	printf("STAT: avg_error = %f p_005 = %f p_095 = %f\n", avg_error, p_005, p_095);
}

double surrogate_estimate(int Gen, double *x, int n)
{
	char netname[256];
	sprintf(netname, "myNN_%03d.txt", Gen);
	double predicted_y = predict_network(netname, x, n);

//	double predicted_y = predict_network("myNN.txt", x, n);

	printf("evaluate_surrogate -> %f\n", predicted_y);

	return predicted_y;
}

double surrogate_error_estimate(int Gen, double *x, int n, double y)
{
	double predicted_y = predict_network("myNN.txt", x, n);

	double ae = fabs(predicted_y-y);
	double re;

	double minv = MIN(fabs(predicted_y),fabs(y));
	double maxv = MAX(fabs(predicted_y),fabs(y));
	if (minv == 0 && maxv == 0)
		re = 0;
	else if (minv == 0)
		re = 100.0*ae / maxv;
	else
		re = 100.0*ae / minv;

//	printf("evaluate_error -> %f vs %f -> %f\n", predicted_y, y, re);

	return re;
}

void destroy_net()
{
//	fann_destroy(ann);
}


