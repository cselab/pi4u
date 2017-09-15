/*
 *  tmcmc_db.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include "engine_tmcmc.h"

void update_curgen_db(cgdb_t *curgen_db, double point[], double F, int dim)
{
	// dim == tmcmc_data.data.Nth;
	int i, pos;

#ifdef USE_TORC
	pthread_mutex_lock(&curgen_db->m);
#endif
	pos = curgen_db->entries;
	curgen_db->entries++;
#ifdef USE_TORC
	pthread_mutex_unlock(&curgen_db->m);
#endif
//	if (curgen_db->entry[pos].point == NULL) curgen_db->entry[pos].point = malloc(tmcmc_data.data.Nth*sizeof(double));
	if (curgen_db->entry[pos].point == NULL) curgen_db->entry[pos].point = malloc(dim*sizeof(double));

	for (i = 0; i < dim; i++) curgen_db->entry[pos].point[i] = point[i];
	curgen_db->entry[pos].F = F;
}

void init_curgen_db(cgdb_t *curgen_db, int popsize)
{
#ifdef USE_TORC
	pthread_mutex_init(&curgen_db->m, NULL);
#endif
	curgen_db->entries = 0;
//	curgen_db->entry = calloc(1, tmcmc_data.data.PopSize*sizeof(cgdbp_t));
	curgen_db->entry = calloc(1, popsize*sizeof(cgdbp_t));
	printf("curgen_db = %p entry = %p\n", curgen_db, curgen_db->entry);
}

void dump_curgen_db(cgdb_t *curgen_db, int Gen, int dim, int t_info[4])
{
//	int PROBDIM = tmcmc_data.data.Nth;
	int pos;
	FILE *fp;
	char fname[256];

#ifdef USE_TORC
//	sprintf(fname, "curgen_db_%d.%d.%d.%d_%03d.txt", t_info[0], t_info[1], t_info[2], t_info[3], Gen);
#else
	sprintf(fname, "curgen_db_%03d_%03d.txt", t_info[0], Gen);
#endif
	fp = fopen(fname, "w");
	for (pos = 0; pos < curgen_db->entries; pos++) {
		if (dim == 2)
			fprintf(fp, "%20.16lf %20.16lf %20.16lf\n",
				curgen_db->entry[pos].point[0], curgen_db->entry[pos].point[1], curgen_db->entry[pos].F);
		else if (dim == 3) 
			fprintf(fp, "%20.16lf %20.16lf %20.16lf %20.16lf\n",
				curgen_db->entry[pos].point[0], curgen_db->entry[pos].point[1], curgen_db->entry[pos].point[2], curgen_db->entry[pos].F);
		else {
			int i;
			
			for (i = 0; i < dim; i++) {
				fprintf(fp, "%20.16lf ", curgen_db->entry[pos].point[i]);
			}
			fprintf(fp, "%20.16lf\n", curgen_db->entry[pos].F);
		}

	}
	fclose(fp);
}
