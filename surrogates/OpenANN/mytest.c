#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>	
#include <math.h>	 
#include <float.h>	// DBL_MAX
#include <string.h>	// strcat

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef MAXLINE
#define MAXLINE 4096
#endif

#ifndef DIM
#define DIM 2
#endif

extern void build_network(char *netname, double *Xin, double *Tin, int D, int F, int N, double *Eout);
extern double predict_network(char *netname, double *Xin, int D);

void construct_nn(char *datafile, int trainsize, int trainoffset, double threshold)
{
	// Create dataset
	const int D = DIM; // number of inputs
	const int F = 1; // number of outputs
	const int N = trainsize; // size of training set
	int real_N;
	int i;

	double *X = calloc(1, N*D*sizeof(double));
	double *T = calloc(1, N*F*sizeof(double));
	double *E = calloc(1, N*sizeof(double));
	int *idx = calloc(1, N*sizeof(double));

	char line[MAXLINE];
	FILE *fp = fopen(datafile, "r");
	if (fp == NULL) {
		printf("Cannot open file!\n"); exit(1);
	}

	for (i = 0; i < trainoffset; i++) {
		fgets(line, MAXLINE, fp);
	}


	int real_i = 0;	//
	for (i = 0; i < N; i++) {
		int j;
		double x[D+F], y;

		fgets(line, MAXLINE, fp);
#if 0
#if DIM == 2
		sscanf(line, "%lf %lf %lf", &x[0], &x[1], &y);
#endif
#else
		int error;
		error = 0;
		char *tok = strtok(line, " ;,\t");
		if (tok == NULL) continue;	// something went wrong
		for (j = 0; j < D+F; j++) {
			if (tok != NULL) {
				x[j] = atof(tok);
			}
			else {
				error = 1;
				break;
			}
			tok = strtok(NULL, " ;,\t");
		}
		if (error) continue;
		y = x[D];
#endif
		if (y < threshold) {
			idx[i] = -1;
			continue;
		} else {
			idx[i] = real_i;
		}

		for (j = 0; j < D; j++) X[real_i*D+j] = x[j];
		T[real_i*F+0] = y;
		real_i++;

	}
	real_N = real_i;

	fclose(fp);

	build_network("myNN.txt", X, T, D, F, real_N, E);

#if VERBOSE
	for (i = 0; i < N; i++) {
		if (idx[i] == -1) 
			printf("%d -> -\n", i);
		else
			printf("%d -> %e\n", i, E[idx[i]]); 
	}
#endif

	free(X);
	free(T);
	free(E);
	free(idx);
}

void evaluate_nn(char *testfile, int testsize, int testoffset, double unused_threshold)
{
	// Create dataset
	const int D = DIM; // number of inputs
	const int F = 1; // number of outputs

	int i, j;

	// Create network

	double mse = 0;
	double mae = 0;
	double mre = 0;
	double sae = 0;
	double sare = 0;

	char line[MAXLINE];
	FILE *fp2 = fopen(testfile, "r");
	for (i = 0; i < testoffset; i++) {
		fgets(line, MAXLINE, fp2);
	}

	#define NBINS   20
	int bins[NBINS+1], bins2[NBINS+1];
	for (i = 0; i < NBINS+1; i++) bins[i]=0;
	for (i = 0; i < NBINS+1; i++) bins2[i]=0;

	int real_testsize = 0;
	for (i = 0; i < testsize; i++) {
		double ae;
		double re;
		double x[D+F], y;

		fgets(line, MAXLINE, fp2);
#if 0
#if DIM==2
		sscanf(line, "%lf %lf %lf", &x[0], &x[1], &y);
#endif
#else
		int error;
		error = 0;
		char *tok = strtok(line, " ;,\t");
		if (tok == NULL) continue;	// something went wrong
		for (j = 0; j < D+F; j++) {
			if (tok != NULL) {
				x[j] = atof(tok);
			}
			else {
				error = 1;
				break;
			}
			tok = strtok(NULL, " ;,\t");
		}
		if (error) continue;
		y = x[D];
#endif

		double predicted_sy = predict_network("myNN.txt", x, D);

		ae = fabs(predicted_sy-y);
		sae += ae;

		double minv = MIN(fabs(predicted_sy),fabs(y));
		double maxv = MAX(fabs(predicted_sy),fabs(y));
		if ((minv == 0)&&(maxv == 0))
			re = 0;
		else if (minv == 0)
			re = 100.0*ae / maxv;
		else
			re = 100.0*ae / minv;

		sare += re;

		double symre = 100.0*ae/(0.5*(fabs(predicted_sy)+fabs(y)));

#if VERBOSE
		int k;
		char x_str[MAXLINE];
		sprintf(x_str, "%f ", x[0]);
		for (k = 1; k < D; k++)
		{
			char x_k[MAXLINE];

			sprintf(x_k, "%f ", x[k]);
			strcat(x_str, x_k);
		}
		printf("out: %s -> %f vs %f %.2f%% %.2f%%\n", x_str, predicted_sy, y, re, symre);
#endif
		if (ae > mae) mae = ae;
		if (re > mre) mre = re;
		mse += ae*ae;

		int cbin = (int)symre/(100.0/NBINS);
		if (cbin > NBINS) cbin = NBINS;
		if (re > 100) cbin = NBINS;
		bins[cbin]++;

		int cbin2 = (int)re/(100.0/NBINS);
		if (cbin2 > NBINS) cbin2 = NBINS;
		if (re > 100) cbin2 = NBINS;
		bins2[cbin2]++;

		real_testsize++;
	}

	mse = mse/(1.0*real_testsize);
	printf("max abs error = %f\n", mae);
	printf("avg abs error = %f\n", sae/(1.0*real_testsize));
	printf("max abs rel error = %.2f%%\n", mre);
	printf("avg abs rel error = %.2f%%\n", sare/(1.0*real_testsize));
	printf("mse = %f\n", mse);

#if 0
	int sbin = 0;
	for (i = 0; i < NBINS+1; i++) sbin+=bins[i];
	printf("\n--SYMRE-\n");
	for (i = 0; i < NBINS; i++)
	{
		printf("[%3.0f - %3.0f) -> %d (%.2f%%)\n", (100.0/NBINS)*i, (100.0/NBINS)*(i+1), bins[i], 100.0*bins[i]/sbin);
	}
	printf("[%3.0f - inf) -> %d (%.2f%%)\n", 100.0, bins[NBINS], 100.0*bins[NBINS]/sbin);
#endif

#if 1
	printf("\n--RE--\n");
	int sbin2 = 0;
	for (i = 0; i < NBINS+1; i++) sbin2+=bins2[i];
	for (i = 0; i < NBINS; i++)
	{
		printf("[%3.0f - %3.0f) -> %d (%.2f%%)\n", (100.0/NBINS)*i, (100.0/NBINS)*(i+1), bins2[i], 100.0*bins2[i]/sbin2);
	}
	printf("[%3.0f - inf) -> %d (%.2f%%)\n", 100.0, bins2[NBINS], 100.0*bins2[NBINS]/sbin2);
#endif
}


int main(int argc, char *argv[])
{
	if (argc < 4)
	{
		printf("usage: [export THRESHOLD=<value>] ./%s <datafile> <trainsize> <trainoffset> [<testfile> <testsize>] [<testoffset>]\n", argv[0]);
		exit(1);
	}

	char *datafile = NULL, *testfile = NULL;
	int trainsize = 0, testsize = 0, trainoffset = 0, testoffset = 0;

	if ((argc == 6)||(argc == 7))
	{
		datafile = argv[1];
		trainsize = atoi(argv[2]);
		trainoffset = atoi(argv[3]);
		testfile = argv[4];
		testsize = atoi(argv[5]);
	}
	else if ((argc == 4)||(argc == 5))
	{
		datafile = argv[1];
		trainsize = atoi(argv[2]);
		trainoffset = atoi(argv[3]);
		testfile = datafile;
		testsize = trainsize;
	}

	if (argc == 7) testoffset = atoi(argv[6]);
	if (argc == 5) testoffset = atoi(argv[4]);

	double threshold = -DBL_MAX;
	char *v = getenv("THRESHOLD");
	if (v != NULL) threshold = atof(v);

	printf("running benchmark for (%d %d) (%d %d) and threshold = %e\n", trainsize, trainoffset, testsize, testoffset, threshold);
	construct_nn(datafile, trainsize, trainoffset, threshold);
	evaluate_nn(testfile, testsize, testoffset, threshold);

	return 0;
}
