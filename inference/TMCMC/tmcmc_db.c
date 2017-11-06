/*
 *  tmcmc_db.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include "engine_tmcmc.h"

void update_full_db(double point[], double F, double *G, int n, int surrogate)
{
	int PROBDIM = data.Nth;
	int i, pos;

	pthread_mutex_lock(&full_db.m);
	pos = full_db.entries;
	full_db.entries++;
	pthread_mutex_unlock(&full_db.m);

	if (full_db.entry[pos].point == NULL) full_db.entry[pos].point = (double *)malloc(data.Nth*sizeof(double));

	for (i = 0; i < PROBDIM; i++) full_db.entry[pos].point[i] = point[i];
	full_db.entry[pos].F = F;
	full_db.entry[pos].nG = n;
	for (i = 0; i < n; i++) full_db.entry[pos].G[i] = G[i];
	full_db.entry[pos].surrogate = surrogate;
}

void init_full_db()
{
	pthread_mutex_init(&full_db.m, NULL);
	full_db.entries = 0;
	full_db.entry = (dbp_t *)calloc(1, data.MaxStages*data.PopSize*sizeof(dbp_t));	
}


void update_curgen_db(double point[], double F, double prior)
{
	int PROBDIM = data.Nth;
	int i, pos;

	pthread_mutex_lock(&curgen_db.m);
	pos = curgen_db.entries;
	curgen_db.entries++;
	pthread_mutex_unlock(&curgen_db.m);

	if (curgen_db.entry[pos].point == NULL) curgen_db.entry[pos].point = (double *)malloc(data.Nth*sizeof(double));

	for (i = 0; i < PROBDIM; i++) curgen_db.entry[pos].point[i] = point[i];
	curgen_db.entry[pos].F = F;
	curgen_db.entry[pos].prior = prior;
}


void init_curgen_db()
{
	pthread_mutex_init(&curgen_db.m, NULL);
	curgen_db.entries = 0;
	curgen_db.entry = (cgdbp_t *)calloc(1, (data.MinChainLength+1)*data.PopSize*sizeof(cgdbp_t));
}


void update_curres_db(double point[EXPERIMENTAL_RESULTS], double F)
{
	int i, pos;
	
#if (EXPERIMENTAL_RESULTS <=0)
	return; 
#endif
	pthread_mutex_lock(&curres_db.m);
	pos = curres_db.entries;
	curres_db.entries++;
	pthread_mutex_unlock(&curres_db.m);

	if (curres_db.entry[pos].point == NULL) curres_db.entry[pos].point = (double *)malloc((EXPERIMENTAL_RESULTS+1)*sizeof(double));

	for (i = 0; i < EXPERIMENTAL_RESULTS; i++) curres_db.entry[pos].point[i] = point[i];
	curres_db.entry[pos].F = F;	
}

void init_curres_db()
{
	pthread_mutex_init(&curres_db.m, NULL);
	curres_db.entries = 0;
	curgen_db.entry = (cgdbp_t *)calloc(1, (data.MinChainLength+1)*data.PopSize*sizeof(cgdbp_t));
}


void print_full_db()
{
	int pos, i;

	printf("=======\n");
	printf("FULL_DB\n");
	for (pos = 0; pos < full_db.entries; pos++) {
		printf("ENTRY %d: POINT(%20.16lf,%20.16lf) F=%20.16lf SG=%d\n",
				pos, full_db.entry[pos].point[0], full_db.entry[pos].point[1],	/* extend it*/
				full_db.entry[pos].F, full_db.entry[pos].surrogate);
		printf("\tG=[");
		for (i = 0; i < full_db.entry[pos].nG-1; i++) printf("%20.16lf,", full_db.entry[pos].G[i]);
		printf("%20.16lf]\n", full_db.entry[pos].G[i]);
	}
	printf("=======\n");
}

void dump_full_db(int Gen)
{
	int PROBDIM = data.Nth;
	int pos;
	FILE *fp;
	char fname[256];

	sprintf(fname, "full_db_%03d.txt", Gen);
	fp = fopen(fname, "w");
	for (pos = 0; pos < full_db.entries; pos++) {
		int i;
		for (i = 0; i < PROBDIM; i++) {
			fprintf(fp, "%20.16lf ", full_db.entry[pos].point[i]);
		}
		fprintf(fp, "%20.16lf\n", full_db.entry[pos].F);
	}
	fclose(fp);
}


void dump_curgen_db(int Gen)
{
	int PROBDIM = data.Nth;
	int pos;
	FILE *fp;
	char fname[256];

	sprintf(fname, "curgen_db_%03d.txt", Gen);
	fp = fopen(fname, "w");
	for (pos = 0; pos < curgen_db.entries; pos++) {
		int i;
			
		for (i = 0; i < PROBDIM; i++) {
			fprintf(fp, "%20.16lf ", curgen_db.entry[pos].point[i]);
		}
		fprintf(fp, "%20.16lf  ", curgen_db.entry[pos].F);
		fprintf(fp, "%20.16lf  ", curgen_db.entry[pos].prior);
		fprintf(fp,"\n");

	}
	fclose(fp);
}

int load_curgen_db(int Gen)
{
	int PROBDIM = data.Nth;	/* peh: be careful here */
	int pos;
	FILE *fp;
	char fname[256];

	sprintf(fname, "curgen_db_%03d.txt", Gen);
	fp = fopen(fname, "r");
	if (fp == NULL) {
		printf("DB file: %s not found!!!\n", fname);
		exit(1); 
		return 1;
	}

	curgen_db.entries = 0;
	char line[1024];
	while (fgets(line, 1024, fp) != NULL)
		curgen_db.entries++;

	fclose(fp);	/* fseek...*/
	fp = fopen(fname, "r");	

	for (pos = 0; pos < curgen_db.entries; pos++) {
		int i;
		for (i = 0; i < PROBDIM; i++) {
			if (curgen_db.entry[pos].point == NULL) curgen_db.entry[pos].point = (double *)malloc(PROBDIM*sizeof(double));
			fscanf(fp, "%lf", &curgen_db.entry[pos].point[i]);
		}
		fscanf(fp, "%lf", &curgen_db.entry[pos].F);
		fscanf(fp, "%lf", &curgen_db.entry[pos].prior);
	}
	fclose(fp);

	return 0;
}

void dump_curres_db(int Gen)
{
	int pos;
	FILE *fp;
	char fname[256];

#if (EXPERIMENTAL_RESULTS <=0)
	return;
#endif
	sprintf(fname, "curres_db_%03d.txt", Gen);
	fp = fopen(fname, "w");
	for (pos = 0; pos < curres_db.entries; pos++) {
		int i;
			
		for (i = 0; i < EXPERIMENTAL_RESULTS; i++) {
			fprintf(fp, "%20.16lf ", curres_db.entry[pos].point[i]);
		}
		fprintf(fp, "%20.16lf\n", curres_db.entry[pos].F);
/*		fprintf(fp, "\n");*/
	}
	fclose(fp);
}
