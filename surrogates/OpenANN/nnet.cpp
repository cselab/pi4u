#include <OpenANN/OpenANN>
#include <OpenANN/io/Logger.h>
#include <OpenANN/io/DirectStorageDataSet.h>
#include <OpenANN/util/Random.h>
#include <Eigen/Core>
#include <iostream>
#include <omp.h>
#include <pthread.h>
using namespace OpenANN;

#ifndef MIN
#define MIN(a,b) ((a)<(b))?(a):(b)
#endif

#ifndef MAX
#define MAX(a,b) ((a)>(b))?(a):(b)
#endif

#ifdef __cplusplus
extern "C" {
#endif

// global data - to be restructured 
#ifndef MAXDIM
#define MAXDIM 64
#endif

// todo: we need to broadcast this data
// todo: thread-safety
struct nnet_data_s
{
	double ymin, ymax, ymed;
	double ymin_sh, ymax_sh;
	double abs_ymax;
	double xmin[MAXDIM], xmax[MAXDIM], xmed[MAXDIM];
	double xmin_sh[MAXDIM], xmax_sh[MAXDIM];
	double amin, amax;
} nnet_data;

static Net *net = NULL;

static char activename[256];
static pthread_mutex_t filemutex = PTHREAD_MUTEX_INITIALIZER;

void save_global(char *name)
{
	int i;

	FILE *fp = fopen(name, "w");
	if (fp == NULL) {
		printf("Could not open file %s for writing!\n", name); exit(1);
	}

	fprintf(fp, "%f\n", nnet_data.ymin);
	fprintf(fp, "%f\n", nnet_data.ymax);
	fprintf(fp, "%f\n", nnet_data.ymed);
	fprintf(fp, "%f\n", nnet_data.ymin_sh);
	fprintf(fp, "%f\n", nnet_data.ymax_sh);
	fprintf(fp, "%f\n", nnet_data.abs_ymax);
	for (i = 0; i < MAXDIM; i++) fprintf(fp, "%f\n", nnet_data.xmin[i]);
	for (i = 0; i < MAXDIM; i++) fprintf(fp, "%f\n", nnet_data.xmax[i]);
	for (i = 0; i < MAXDIM; i++) fprintf(fp, "%f\n", nnet_data.xmed[i]);
	for (i = 0; i < MAXDIM; i++) fprintf(fp, "%f\n", nnet_data.xmin_sh[i]);
	for (i = 0; i < MAXDIM; i++) fprintf(fp, "%f\n", nnet_data.xmax_sh[i]);
	fprintf(fp, "%f\n", nnet_data.amin);
	fprintf(fp, "%f\n", nnet_data.amax);	
	fclose(fp);

	strcpy(activename, name);
}

void load_global(char *name)
{
	int i;
	if (strcmp(name, activename) == 0) return;

	if (pthread_mutex_trylock(&filemutex) == EBUSY)	// this is tricky and/but I do not like it much.
	{
		pthread_mutex_lock(&filemutex);
		pthread_mutex_unlock(&filemutex);
		return;
	}

	FILE *fp = fopen(name, "r");
	if (fp == NULL) {
		printf("Could not open file %s for reading!\n", name); exit(1);
	}
	fscanf(fp, "%lf\n", &nnet_data.ymin);
	fscanf(fp, "%lf\n", &nnet_data.ymax);
	fscanf(fp, "%lf\n", &nnet_data.ymed);
	fscanf(fp, "%lf\n", &nnet_data.ymin_sh);
	fscanf(fp, "%lf\n", &nnet_data.ymax_sh);
	fscanf(fp, "%lf\n", &nnet_data.abs_ymax);
	for (i = 0; i < MAXDIM; i++) fscanf(fp, "%lf\n", &nnet_data.xmin[i]);
	for (i = 0; i < MAXDIM; i++) fscanf(fp, "%lf\n", &nnet_data.xmax[i]);
	for (i = 0; i < MAXDIM; i++) fscanf(fp, "%lf\n", &nnet_data.xmed[i]);
	for (i = 0; i < MAXDIM; i++) fscanf(fp, "%lf\n", &nnet_data.xmin_sh[i]);
	for (i = 0; i < MAXDIM; i++) fscanf(fp, "%lf\n", &nnet_data.xmax_sh[i]);
	fscanf(fp, "%lf\n", &nnet_data.amin);
	fscanf(fp, "%lf\n", &nnet_data.amax);	
	fclose(fp);

	strcpy(activename, name);
	pthread_mutex_unlock(&filemutex);
}

// the following options should be included in a configuration file
#define FULLY_CONNECTED
#define HIDDEN_NEURONS	32

#define SCALE_OUTPUT
#define SCALE_INPUT

double mapminmax_apply(double x, double min_v, double max_v, double shift_v, double min_a, double max_a)
{
#ifndef SCALE_OUTPUT
	return x;
#else
	x = x - shift_v;
	x = (max_a-min_a)*(x-min_v)/(max_v-min_v) + (min_a);
	return x;
#endif
}

double mapminmax_reverse(double x, double min_v, double max_v, double shift_v, double min_a, double max_a)
{
#ifndef SCALE_OUTPUT
	return x;
#else
	x = (max_a-min_a)*(x-min_v)/(max_v-min_v) + (min_a);
	x = x + shift_v;
	return x;
#endif
}

// netname: file to store NN
// Xin: evaluation points
// Tin: target values
// D: number of inputs
// F: number of outputs (1)
// N: size of training set
// Eout: error for each point
void build_network(char *netname, double *Xin, double *Tin, int D, int F, int N, double *Eout)
{
	if (D > MAXDIM) {
		printf("Number of inputs (D) bigger than MAXDIM. Modify and recompile!\n");
		exit(1);
	}

	// Create dataset
	Eigen::MatrixXd X(N, D); // inputs
	Eigen::MatrixXd T(N, F); // desired outputs (targets)

	// compute min-max values
	nnet_data.ymin = nnet_data.ymax = Tin[0];
	for (int j = 0; j < D; j++) {
		nnet_data.xmin[j] = nnet_data.xmax[j] = Xin[0*D+j];
	}

	for (int i = 0; i < N; i++) {
		double xj, y;

		for (int j = 0; j < D; j++) {
			xj = Xin[i*D+j];
			if (xj < nnet_data.xmin[j]) nnet_data.xmin[j] = xj;
			if (xj > nnet_data.xmax[j]) nnet_data.xmax[j] = xj;
		}

		y = Tin[i];
		if (y < nnet_data.ymin) nnet_data.ymin = y;
		if (y > nnet_data.ymax) nnet_data.ymax = y;
		if (fabs(y) > nnet_data.abs_ymax) nnet_data.abs_ymax = fabs(y);

	}

	nnet_data.amin = -1;	// fixed
	nnet_data.amax = +1;	// fixed
	nnet_data.ymed = 0.5*(nnet_data.ymax + nnet_data.ymin);
	nnet_data.ymin_sh = nnet_data.ymin - nnet_data.ymed;
	nnet_data.ymax_sh = nnet_data.ymax - nnet_data.ymed;
#if VERBOSE
	printf("y = [%lf %lf] -> [%lf %lf] ymed = %lf\n", ymin, ymax, ymin_sh, ymax_sh, ymed);
#endif

	for (int j = 0; j < D; j++) {
		nnet_data.xmed[j] = 0.5*(nnet_data.xmax[j]+nnet_data.xmin[j]);
		nnet_data.xmin_sh[j] = nnet_data.xmin[j] - nnet_data.xmed[j]; 
		nnet_data.xmax_sh[j] = nnet_data.xmax[j] - nnet_data.xmed[j]; 
#if VERBOSE
		printf("%d: x = [%lf %lf] -> [%lf %lf] xmed = %lf\n", j, nnet_data.xmin[j], nnet_data.xmax[j], nnet_data.xmin_sh[j], nnet_data.xmax_sh[j], nnet_data.xmed[j]);
#endif
	}


	// prepare dataset
	for (int i = 0; i < N; i++) {
		double x[D], y;

		for (int j = 0; j < D; j++) x[j] = Xin[i*D+j];
#ifdef SCALE_INPUT
		for (int j = 0; j < D; j++) {
			x[j] = mapminmax_apply(x[j], nnet_data.xmin_sh[j], nnet_data.xmax_sh[j], nnet_data.xmed[j], nnet_data.amin, nnet_data.amax);
#if VERBOSE
			printf("x %d/%d: from %f to %f\n", i, j, Xin[i*D+j], x[j]);
#endif
		}
#endif
		for (int j = 0; j < D; j++) X(i,j) = x[j];

		y = Tin[i];
		y = mapminmax_apply(y, nnet_data.ymin_sh, nnet_data.ymax_sh, nnet_data.ymed, nnet_data.amin, nnet_data.amax);
		T(i,0) = y;
#if VERBOSE
		printf("y %d: from %f to %f\n", i, Tin[i], y);
#endif
	}

	// prepare training data for OpenANN
	DirectStorageDataSet trainingSet(&X, &T);

	// Make the result repeatable
	RandomNumberGenerator().seed(0);

	// Create network
	net = new Net;

	// input +  hidden + output 
	net->inputLayer(D);

	// todo: read NN configuration from input file 
#ifdef FULLY_CONNECTED
	net->fullyConnectedLayer(HIDDEN_NEURONS, LOGISTIC);
#else
	net->extremeLayer(512, LOGISTIC);
#endif
	net->outputLayer(F, LINEAR);

	// Add training set
	net->trainingSet(trainingSet);

	// Set stopping conditions
	StoppingCriteria stop;
#ifdef FULLY_CONNECTED
	stop.minimalValueDifferences = 0;
#else
	stop.minimalValueDifferences = 1e-10;
#endif
	stop.maximalIterations = 1000;

	// Train network, i.e. minimize sum of squared errors (SSE) with
	// Levenberg-Marquardt optimization algorithm until the stopping criteria
	// are satisfied.
	double t0 = omp_get_wtime();
	train(*net, "LMA", MSE, stop);
	double t1 = omp_get_wtime();
	printf("training time = %f secs\n", t1-t0);

	{
	#define NBINS   20
	int bins[NBINS+1];
	for (int i = 0; i < NBINS+1; i++) bins[i]=0;

	for (int i = 0; i < N; i++) {
		double ae, re, symre;

		Eigen::VectorXd output = (*net)(trainingSet.getInstance(i));
		double predicted_y = output(0);

		double y = Tin[i];
		double predicted_sy = mapminmax_reverse(predicted_y, nnet_data.amin, nnet_data.amax, nnet_data.ymed, nnet_data.ymin_sh, nnet_data.ymax_sh);

		ae = fabs(predicted_sy-y);

		double minv = MIN(fabs(predicted_sy),fabs(y));
		double maxv = MAX(fabs(predicted_sy),fabs(y));
		if (minv == 0 && maxv == 0)
			re = 0;
		else if (minv == 0)
			re = 100.0*ae / maxv;
		else
			re = 100.0*ae / minv;

		symre = 100.0*ae/(0.5*(fabs(predicted_sy)+fabs(y)));

#if VERBOSE
//		double x[D], y;
//		for (int j = 0; j < D; j++) x[j] = X(i,j);

		printf("%d -> %f vs %f %.2f%% %.2f%%\n", i, predicted_sy, y, re, symre);
#endif

		double error;	// re or symre
		if (1)
			error = re;
		else
			error = symre;

		if (Eout != NULL) Eout[i] = error;		

		int cbin = (int)error/(100.0/NBINS);
		if (cbin > NBINS) cbin = NBINS;
		if (error > 100.0) cbin = NBINS;
		bins[cbin]++;
	}

	// this is optional, but we keep it
	int sbin = 0;
	for (int i = 0; i < NBINS+1; i++) sbin+=bins[i];
	printf("\n- SELF-EVALUATION: HISTOGRAM -\n");
	for (int i = 0; i < NBINS; i++)
	{
		printf("[%3.0f - %3.0f) -> %d (%.2f%%)\n", (100.0/NBINS)*i, (100.0/NBINS)*(i+1), bins[i], 100.0*bins[i]/sbin);
	}
	printf("[%3.0f - inf) -> %d (%.2f%%)\n", 100.0, bins[NBINS], 100.0*bins[NBINS]/sbin);
	}

	net->save(netname);
#if 1
	char global_name[256];
	sprintf(global_name, "scale_%s", netname);
	save_global(global_name);
#endif
	delete net;
	net = NULL;
}

// netname: file with stored NN
// Xin: evaluation point
// D: number of inputs
double predict_network(char *netname, double *Xin, int D)
{
	// Create network
	Net *lnet = new Net;

	// this should be loaded once
	lnet->load(netname);
#if 1
	char global_name[256];
	sprintf(global_name, "scale_%s", netname);
	load_global(global_name);
#endif

	double x[D];
	for (int j = 0; j < D; j++) x[j] = Xin[j];

#ifdef SCALE_INPUT
	for (int j = 0; j < D; j++) x[j] = mapminmax_apply(x[j], nnet_data.xmin_sh[j], nnet_data.xmax_sh[j], nnet_data.xmed[j], nnet_data.amin, nnet_data.amax);
#endif

	Eigen::VectorXd Xt(D);
	for (int j = 0; j < D; j++) Xt(j) = x[j];

	Eigen::VectorXd output = (*lnet)(Xt);
	double predicted_y = output(0);
	double predicted_sy = mapminmax_reverse(predicted_y, nnet_data.amin, nnet_data.amax, nnet_data.ymed, nnet_data.ymin_sh, nnet_data.ymax_sh);

	delete lnet;
	return predicted_sy;
}

#ifdef __cplusplus
}
#endif

