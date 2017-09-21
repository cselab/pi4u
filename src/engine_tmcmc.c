/*
 *  engine_tmcmc.c
 *  Pi4U 
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <assert.h>
#include "engine_tmcmc.h"
#include <dlfcn.h>
//#define _STEALING_

#ifndef USE_TORC
#include <sys/time.h>
static double torc_gettime()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}
#endif

void data_init(tmcmc_data_t *tmcmc_data, int Nth, int MaxStages, int PopSize, double *lb, double *ub)
{
	int i;

	/* DATA: user's input parameters */
	/* instead of read_data() */

	tmcmc_data->data.Nth = Nth;
	tmcmc_data->data.MaxStages = MaxStages;
	tmcmc_data->data.PopSize = PopSize;

	tmcmc_data->data.lb = -6;	// Default LB, same for all
	tmcmc_data->data.ub = +6;	// Default UB, same for all
 
	tmcmc_data->data.lowerbound = malloc(tmcmc_data->data.Nth*sizeof(double));
	tmcmc_data->data.upperbound = malloc(tmcmc_data->data.Nth*sizeof(double));

	for (i = 0; i < tmcmc_data->data.Nth; i++) {
		tmcmc_data->data.lowerbound[i] = lb[i];
		tmcmc_data->data.upperbound[i] = ub[i];
	}

	tmcmc_data->data.TolCOV = 1.0;	// 0.25, 0.5
	tmcmc_data->data.bbeta = 0.04;
	tmcmc_data->data.seed = 280675;

	tmcmc_data->data.options.MaxIter = 100;	// 
	tmcmc_data->data.options.Tol = 1e-6;
	tmcmc_data->data.options.Display = 0;

	tmcmc_data->data.iplot = 0;	// gnuplot
	tmcmc_data->data.idump = 1;	// dump db

	tmcmc_data->data.Num = malloc(tmcmc_data->data.MaxStages*sizeof(double));
	for (i = 0; i < tmcmc_data->data.MaxStages; i++) {
		tmcmc_data->data.Num[i] = tmcmc_data->data.PopSize; // default DATANUM
	}
	tmcmc_data->data.LastNum = tmcmc_data->data.PopSize; // DATANUM;

#if 1
	init_curgen_db(&tmcmc_data->curgen_db, tmcmc_data->data.PopSize);
#endif

	/* RUNINFO: running state */
	tmcmc_data->runinfo.CoefVar = calloc(1, tmcmc_data->data.MaxStages*sizeof(double));
	tmcmc_data->runinfo.p = calloc(1, tmcmc_data->data.MaxStages*sizeof(double));
	tmcmc_data->runinfo.currentuniques = calloc(1, tmcmc_data->data.MaxStages*sizeof(double));
	tmcmc_data->runinfo.logselection = calloc(1, tmcmc_data->data.MaxStages*sizeof(double));
	tmcmc_data->runinfo.acceptance = calloc(1, tmcmc_data->data.MaxStages*sizeof(double));

	double *SSmem = (double *)calloc(1, tmcmc_data->data.Nth*tmcmc_data->data.Nth*sizeof(double));
	tmcmc_data->runinfo.SS = (double **)malloc(tmcmc_data->data.Nth*sizeof(double *));
	for (i = 0; i < tmcmc_data->data.Nth; i++) {
		tmcmc_data->runinfo.SS[i] = SSmem + i*tmcmc_data->data.Nth; //&SSmem[i*data.Nth];
	}

	tmcmc_data->runinfo.meantheta = calloc(1, tmcmc_data->data.MaxStages*sizeof(double *));
	for (i = 0; i < tmcmc_data->data.MaxStages; i++) {
		tmcmc_data->runinfo.meantheta[i] = calloc(1, tmcmc_data->data.Nth*sizeof(double));
	}

	tmcmc_data->runinfo.Gen = 0;
	tmcmc_data->runinfo.CoefVar[0] = 10;

	/* already zero */
	tmcmc_data->runinfo.p[0] = 0;
	for (i = 0; i < tmcmc_data->data.MaxStages; i++) {
		tmcmc_data->runinfo.currentuniques[i] = 0;
		tmcmc_data->runinfo.logselection[i] = 0.0;
		tmcmc_data->runinfo.acceptance[i] = 0.0;
	}
}

void check_for_exit(tmcmc_data_t *tmcmc_data)
{
	int val, exitgen = -1;
	char *s;

	s = (char *) getenv("EXIT_GEN");
	if (s != 0 && sscanf(s, "%d", &val) == 1 && val >= 0)
		exitgen = val;

	if (exitgen == tmcmc_data->runinfo.Gen) {
		printf("Read Exit Envrironment Variable!!!\n");
#ifdef USE_TORC
		torc_finalize();
#endif
		exit(1);
	}

	FILE *fp;
	fp = fopen("exit.txt", "r");
	if (fp != NULL) {
		printf("Found Exit File!!!\n");
		unlink("exit.txt");
#ifdef USE_TORC
		torc_finalize();
#endif
		exit(1);
	}
}

#if 1	/* TORC-BASED DATA MANAGEMENT */
void torc_update_curgen_db_task(double point[], double *pF, tmcmc_data_t *tmcmc_data, int *tmcmc_owner)
{
	double F = *pF;
//	printf("%d torc_update_curgen_db_task: %f %f %f %p %d-%d\n", me, point[0], point[1], F, tmcmc_data, tmcmc_owner[0], tmcmc_owner[1]); fflush(0);

	update_curgen_db(&tmcmc_data->curgen_db, point, F, tmcmc_data->data.Nth);
}

void torc_update_curgen_db(double point[], double F, int dim, tmcmc_data_t *tmcmc_data, int *tmcmc_owner)
{
#ifdef USE_TORC
	int me = torc_node_id();
#else
	int me = 0;
#endif
	//if (me == 0) {
	//printf("%d torc_update_curgen_db: %f %f %f %p %d-%d\n", me, point[0], point[1], F, tmcmc_data, tmcmc_owner[0], tmcmc_owner[1]); fflush(0);
	//}

	if (me == tmcmc_owner[0]) {	
		update_curgen_db(&tmcmc_data->curgen_db, point, F, dim);
		return;
	}

#ifdef USE_TORC
	// xxx: add some reference to tmcmc_data, 0 -> tmcmc owner 
	torc_create_direct(tmcmc_owner[1], torc_update_curgen_db_task, 4,		/* message to the database manager (separate process?) or direct execution by server thread */
		dim, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_LONG, CALL_BY_VAD,
		2, MPI_INT, CALL_BY_COP,
		point, &F, tmcmc_data, tmcmc_owner);
	torc_waitall3();	// wait without releasing the worker ...
#endif
}
#endif


/*** TASK MANAGEMENT ***/

#if 0
typedef double (*feval_t)(double *x, long n, void *output, int *info);
feval_t fitfun;
#else
double fitfun(double *TP, long N, void *output, int *info);
#endif

void taskfun(double /*const*/ *x, int *pN, double *res, int winfo[4])
{
	double f;
	int N = *pN;
//	printf("taskfun {%d,%d,%d,%d}\n", winfo[0], winfo[1], winfo[2], winfo[3]);

	inc_nfc();	// increment function call counter

	f = fitfun(x, N, (void *)NULL, winfo);
	*res = f;
	return;
}

double F(double *TP, int *pn)	/* for PNDL */
{
	double gres;

	taskfun(TP, pn, &gres, NULL);
//	gres = posterior(TP, PROBDIM, gres);

	return gres;
}

void evaluate_F(double point[], int dim, double *Fval, int gen_id, int chain_id, int step_id, int ntasks)
{
	int i;
	double G[64], F;	// maxtasks
	int winfo[4];

	for (i = 0; i < ntasks; i++) {
		winfo[0] = gen_id;
		winfo[1] = chain_id;
		winfo[2] = step_id;
		winfo[3] = i;

		if (ntasks == 1) {
			taskfun(point, &dim, &G[i], winfo);
		}
		else {
			exit(1);
#if 0
			torc_create(-1, taskfun, 4,
				dim, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_RES,
				4, MPI_INT, CALL_BY_COP,
				point, &dim, &G[i], winfo);
#endif
		}

	}
#ifdef USE_TORC
	if (ntasks > 1)
		torc_waitall();
#endif

	/* compute F from G's - how? : F fitness function = distance from ground truth */
	F = 0.0;
	for (i = 0; i < ntasks; i++) {
		F += G[i];
	}
	F = F/ntasks;

	//update_full_db(point, F, G, ntasks, 0); // not surrogate 
	//...torc_update_full_db(point, F, G, ntasks, 0); // not surrogate 

	*Fval = F;
}

void initchaintask(double in_tparam[], int *pdim, double *out_tparam, int winfo[4], tmcmc_data_t *tmcmc_data, int *tmcmc_owner)
{
	int i;
	int dim = *pdim;
//	int tmcmc_owner = *ptmcmc_owner;
	int gen_id, chain_id;
	gen_id = winfo[0];
	chain_id = winfo[1];
	
	double point[dim], fpoint;

	for (i = 0; i < dim; i++)
		point[i] = in_tparam[i];

	evaluate_F(point, dim, &fpoint, gen_id, chain_id, 0, 1);	// 12 for namd

	/* update current db entry */
	//if (torc_node_id() > 0) {
	//	printf("XXX tmcmc_data = %p\n", tmcmc_data);
	//}
	torc_update_curgen_db(point, fpoint, dim, tmcmc_data, tmcmc_owner);
	*out_tparam = fpoint;	// currently not required, the result is already in the db

	return;
}

void compute_candidate(double candidate[], double leader[], int dim, double var, double *bSS, double *lowerbound, double *upperbound)
{
//	tmcmc_data_t *tmcmc_data = 0;	// xxx
//	int i, j;
	int i;

//	double bSS[dim*dim];
//	for (i = 0; i < dim; i++)
//		for (j = 0; j < dim; j++)
//			printf("bSS[%d][%d]=%lf\n", i, j, bSS[i*dim+j]);

retry:
	mvnrnd(leader, (double *)bSS, candidate, dim);
	for (i = 0; i < dim; i++) {
		if (isnan(candidate[i])) {
			printf("!!!!  isnan in candidate point!\n");
			exit(1);
			break;
		}
		//if ((candidate[i] < tmcmc_data->data.lowerbound[i])||(candidate[i] > tmcmc_data->data.upperbound[i])) break;
		if ((candidate[i] < lowerbound[i])||(candidate[i] > upperbound[i])) break;
	}
	if (i < dim) goto retry;
}

void chaintask(double in_tparam[], int *pdim, int *pnsteps, double *out_tparam, int winfo[4], double *ppj, double *bSS, double *lowerbound, double *upperbound,
tmcmc_data_t *tmcmc_data, int *tmcmc_owner)
{
	int i,step;
	int dim = *pdim;
//	int tmcmc_owner = *ptmcmc_owner;
	int nsteps = *pnsteps;
	int gen_id = winfo[0];
	int chain_id = winfo[1];
	
	double leader[dim], fleader, fpc_leader;		// fold
	double candidate[dim], fcandidate, fpc_candidate;	// fnew

	for (i = 0; i < dim; i++) leader[i] = in_tparam[i]; //chainwork->in_tparam[i];	// get initial leader
	fleader = *out_tparam; //chainwork->out_tparam[0];					// and its value
	fpc_leader = posterior(leader, dim, fleader);
//	double pj = tmcmc_data.runinfo.p[tmcmc_data.runinfo.Gen];
	double pj = *ppj;

	for (step = 0; step < nsteps; step++) {

		compute_candidate(candidate, leader, dim, 1, bSS, lowerbound, upperbound); //bbeta*SS);	// multivariate gaussian(center, var) for each direction

		/* evaluate fcandidate (NAMD: 12 points) */
		evaluate_F(candidate, dim, &fcandidate, gen_id, chain_id, step, 1);	// this can spawn many tasks
		fpc_candidate = posterior(candidate, dim, fcandidate);

		/* Decide */
		double prior_candidate = priorpdf(candidate, dim, lowerbound, upperbound);	// from PanosA 
		double prior_leader = priorpdf(leader, dim, lowerbound, upperbound);
		double L = exp((prior_candidate-prior_leader)+(fpc_candidate-fpc_leader)*pj);
		if (L > 1) L = 1;
		double P = uniformrand(0,1);
		if (P < L) {
			for (i = 0; i < dim; i++) leader[i] = candidate[i];	// new leader! 
			fleader = fcandidate;
			fpc_leader = fpc_candidate;
			torc_update_curgen_db(leader, fleader, dim, tmcmc_data, tmcmc_owner);
		}
		else {
			//increase counter or add the leader again in curgen_db
			torc_update_curgen_db(leader, fleader, dim, tmcmc_data, tmcmc_owner);
		}
	}

	return;
}

#if 1
typedef struct {
	int idx;
	int nsel;
	double F;
} sort_t;

int compar_desc(const void* p1, const void* p2)
{
	int dir = +1;   // -1: ascending order, +1: descending order
	sort_t *s1 = (sort_t *) p1;
	sort_t *s2 = (sort_t *) p2;

	if (s1->nsel < s2->nsel) return dir;
	if (s1->nsel > s2->nsel) return -dir;
//	if (s1->nsel == s2->nsel) return 0;
	return 0;
}
#endif


int prepare_newgen(tmcmc_data_t *tmcmc_data, int nchains, cgdbp_t *leaders)
{
	/* process curgen_db -> calculate statitics */
	/* compute probs based on F values */
	/* draw new samples (nchains or user-specified) */
	/* find unique samples: fill the (new) leaders table */
	/* count how many times they appear -> nsteps */
	/* return the new sample size (number of chains) */

	int i, p;
	int newchains; // = nchains;
	int dim = tmcmc_data->data.Nth;

	int n = tmcmc_data->curgen_db.entries;

	double fj[n];
	unsigned int sel[n];

	double **g_x;
	g_x = (double **)malloc(dim*sizeof(double *));
	for (i = 0; i < dim; i++)
		g_x[i] = (double *)malloc(n*sizeof(double));

	{//start block
	double **x = g_x;

	for (p = 0; p < dim; p++) {
		for (i = 0; i < n; i++) {
			x[p][i] = tmcmc_data->curgen_db.entry[i].point[p];
		}
	}

	double meanx[dim], stdx[dim];
	for (p = 0; p < dim; p++) {
		meanx[p] = compute_mean(x[p], n);
		stdx[p] = compute_std(x[p], n, meanx[p]);
	}

	printf("CURGEN DB (COMPLE) %d\n", tmcmc_data->runinfo.Gen);
	print_matrix("means", meanx, dim);
	print_matrix("std", stdx, dim);
	}//end block

	{//start block
	double **x = g_x;
	int un = 0, unflag, j;

	for (p = 0; p < dim; p++) {
		x[p][un] = tmcmc_data->curgen_db.entry[0].point[p];	// un==0
	}
	un++;
	for (i = 1; i < n; i++) {
		double xi[dim];
		for (p = 0; p < dim; p++) {
			xi[p] = tmcmc_data->curgen_db.entry[i].point[p];
		}
		unflag = 1;	// is this point unique?
		for (j = 0; j < un; j++) {	// check all the previous unique points
			int compflag;
			compflag = 1;	// 
			for (p = 0; p < dim; p++) {
				if (fabs(xi[p]-x[p][j]) > 1e-6) {
				//if (xi[p] != x[p][j]) {
					compflag = 0;	// they differ
					break;
				}
			}
			
			if (compflag == 1) { 
				unflag = 0;	// not unique, just found it in the unique points table
				break;
			}
		}
		if (unflag) {	// unique, put it in the table
			for (p = 0; p < dim; p++) {
				x[p][un] = xi[p];
			}
			un++;
		}
	} // end block

	tmcmc_data->runinfo.currentuniques[tmcmc_data->runinfo.Gen] = un; //+ 1;
	tmcmc_data->runinfo.acceptance[tmcmc_data->runinfo.Gen] = (1.0*tmcmc_data->runinfo.currentuniques[tmcmc_data->runinfo.Gen])/tmcmc_data->data.Num[tmcmc_data->runinfo.Gen]; // check this

	double meanx[dim], stdx[dim];
	for (p = 0; p < dim; p++) {
		meanx[p] = compute_mean(x[p], un);
		stdx[p] = compute_std(x[p], un, meanx[p]);
	}

	printf("CURGEN DB (UNIQUE) %d: [un = %d]\n", tmcmc_data->runinfo.Gen, un); // + 1);
	print_matrix("means", meanx, dim);
	print_matrix("std", stdx, dim);
	} // end block

	for (i = 0; i < n; i++) fj[i] = tmcmc_data->curgen_db.entry[i].F;	// separate point from F ?
	calculate_statistics(tmcmc_data, fj, n, tmcmc_data->data.Num[tmcmc_data->runinfo.Gen], tmcmc_data->runinfo.Gen, sel);

	newchains = 0;
	for (i = 0; i < n; i++) {
		if (sel[i] != 0) newchains++;
	}

	sort_t list[n];
	for (i = 0; i < n; i++) {
		list[i].idx = i;
		list[i].nsel = sel[i];
		list[i].F = tmcmc_data->curgen_db.entry[i].F;
	}

#if VERBOSE
	printf("Points before\n");
	for (i = 0; i < n; i++) {
		printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
	}
#endif
	qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
	printf("Points after\n");
	for (i = 0; i < n; i++) {
		printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
	}
#endif
	int ldi;	// leader index
	ldi = 0;
	for (i = 0; i < n; i++) {	// newleader
		if (list[i].nsel != 0) {
			int idx = list[i].idx;
			for (p = 0; p < dim; p++) {
				leaders[ldi].point[p] = tmcmc_data->curgen_db.entry[idx].point[p];
			}
			leaders[ldi].F = tmcmc_data->curgen_db.entry[idx].F;
			leaders[ldi].nsel = list[i].nsel;
			ldi++;
		}
	}

	for (i = 0; i < newchains; i++) leaders[i].queue = -1;	// rr

#if VERBOSE
	printf("Leaders before\n");
	for (i = 0; i < newchains; i++) {
		printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
	}
#endif

	/* cool and greedy partitioning ala Panos-- ;-) */

#ifdef USE_TORC
	int nworkers = torc_num_workers();
#else
	int nworkers = 1;
#endif
	int *workload = calloc(1, nworkers*sizeof(int));	// workload[1..workers] = 0

	for (i = 0; i < newchains; i++) {
		int least_loader_worker = compute_min_idx_i(workload, nworkers);
		leaders[i].queue = least_loader_worker;
		workload[least_loader_worker] += leaders[i].nsel;
	}

	print_matrix_i("initial workload", workload, nworkers);
	free(workload);

#if VERBOSE
	printf("Leaders after\n");
	for (i = 0; i < newchains; i++) {
		printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
	}
#endif

	{//start block
//	double x[dim][n];
	double **x = g_x;
	for (i = 0; i < newchains; i++) {
		for (p = 0; p < dim; p++) {
			x[p][i] = leaders[i].point[p];
		}
	}
	
	double meanx[dim], stdx[dim];
	for (p = 0; p < dim; p++) {
		meanx[p] = compute_mean(x[p], newchains);
		stdx[p] = compute_std(x[p], newchains, meanx[p]);
	}

	printf("CURGEN DB (LEADER) %d: [nlead=%d]\n", tmcmc_data->runinfo.Gen, newchains);
	print_matrix("means", meanx, dim);
	print_matrix("std", stdx, dim);
	}//end block

	tmcmc_data->curgen_db.entries = 0;	// reset curgen db
	printf("calculate_statistics: newchains=%d\n", newchains);

	for (i = 0; i < dim; i++) free(g_x[i]);
	free(g_x);

	return newchains;
}


void call_gsl_rand_init(int *pseed)
{
	int seed = *pseed; 
//	printf("CALLING gsl_rand_init(%d) on node %d\n", tmcmc_data.data.seed, torc_node_id()); fflush(0);
//	gsl_rand_init(tmcmc_data.data.seed);
//	printf("CALLING gsl_rand_init(%d) on node %d\n", seed, torc_node_id()); fflush(0);
	gsl_rand_init(seed);
	printf("after\n"); fflush(0);
}

void spmd_gsl_rand_init(int seed)
{
#ifdef USE_TORC
	int i;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, call_gsl_rand_init, 1,
					1, MPI_INT, CALL_BY_COP, &seed);
	}
	torc_waitall();
#else
	call_gsl_rand_init(&seed);
#endif
}

void call_print_matrix_2d(tmcmc_data_t *tmcmc_data)
{
	printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	print_matrix_2d("runinfo.SS", tmcmc_data->runinfo.SS, tmcmc_data->data.Nth, tmcmc_data->data.Nth);
	printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
}

//void spmd_print_matrix_2d()
//{
//	int i;
//	for (i = 0; i < torc_num_nodes(); i++) {
//		torc_create_ex(i*torc_i_num_workers(), 1, call_print_matrix_2d, 0);
//	}
//	torc_waitall();
//}

void call_update_gdata()	// step for p[j]
{
//	MPI_Bcast(tmcmc_data.runinfo.SS[0], tmcmc_data.data.Nth*tmcmc_data.data.Nth, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Bcast(tmcmc_data.runinfo.p, tmcmc_data.data.MaxStages, MPI_DOUBLE, 0, MPI_COMM_WORLD);	// just p[Gen]
//	MPI_Bcast(&tmcmc_data.runinfo.Gen, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void spmd_update_gdata()	// step
{
#ifdef USE_TORC
	int i;
	if (torc_num_nodes() == 1) return;
	for (i = 0; i < torc_num_nodes(); i++) {
		torc_create_ex(i*torc_i_num_workers(), 1, call_update_gdata, 0);
	}
	torc_waitall();
#endif
}

void tmcmc(double *res, int tmcmc_info[4], int Nth, int MaxStages, int PopSize, double *lb, double *ub)
{
	tmcmc_data_t tmcmc_data_s;
	tmcmc_data_t *tmcmc_data = &tmcmc_data_s;

	int i;
	double t0, gt0, gt1;
	int winfo[4];
	int nchains; // was below
	int tmcmc_owner[2];	// node and worker
#ifdef USE_TORC
	tmcmc_owner[0] = torc_node_id();
	tmcmc_owner[1] = torc_worker_id();
#else
	tmcmc_owner[0] = 0;
	tmcmc_owner[1] = 0;
#endif

	printf("tmcmc: tmcmc_data = %p on %d-%d\n", tmcmc_data, tmcmc_owner[0], tmcmc_owner[1]);  fflush(0);

//	data_init(tmcmc_data);
	data_init(tmcmc_data, Nth, MaxStages, PopSize, lb, ub);


	spmd_gsl_rand_init(tmcmc_data->data.seed);

	tmcmc_data->curgen_db.entries = 0; // peh+
	gt0 = t0 = torc_gettime();

	nchains = tmcmc_data->data.Num[0];
	double out_tparam[tmcmc_data->data.PopSize];	// nchains
	for (i = 0; i < nchains; i++) {
		winfo[0] = tmcmc_data->runinfo.Gen;
		winfo[1] = i;
		winfo[2] = -1;
		winfo[3] = -1;

		int d;
		double in_tparam[tmcmc_data->data.Nth];
		for (d = 0; d < tmcmc_data->data.Nth; d++) {
			in_tparam[d] = uniformrand(0,1);
			in_tparam[d] *= (tmcmc_data->data.upperbound[d]-tmcmc_data->data.lowerbound[d]);
			in_tparam[d] += tmcmc_data->data.lowerbound[d];
		}

#ifdef USE_TORC
		torc_create(-1, initchaintask, 6,
			tmcmc_data->data.Nth, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			4, MPI_INT, CALL_BY_COP,
			1, MPI_LONG, CALL_BY_VAD,
			2, MPI_INT, CALL_BY_COP,
			in_tparam, &tmcmc_data->data.Nth, &out_tparam[i], winfo, tmcmc_data, tmcmc_owner);
#else
		initchaintask(in_tparam, &tmcmc_data->data.Nth, &out_tparam[i], winfo, tmcmc_data, tmcmc_owner);
#endif

	}

#ifdef USE_TORC
#ifdef _STEALING_
	torc_enable_stealing();
#endif
	torc_waitall();
#ifdef _STEALING_
	torc_disable_stealing();
#endif
#endif

	gt1 = torc_gettime();
	int g_nfeval = get_nfc();
	printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", tmcmc_data->runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
	reset_nfc();

	{
	*res = 1.0;
	//return;
	}

//	print_full_db();
//	print_curgen_db();
	if (tmcmc_data->data.idump) dump_curgen_db(&tmcmc_data->curgen_db, tmcmc_data->runinfo.Gen, tmcmc_data->data.Nth, tmcmc_info);

	// save here
	//static
	cgdbp_t *leaders; //[MAXCHAINS];
	leaders = calloc(1, tmcmc_data->data.PopSize*sizeof(cgdbp_t));
	for (i = 0; i < tmcmc_data->data.PopSize; i++) {
		leaders[i].point = calloc(1, tmcmc_data->data.Nth*sizeof(double));
	}

	nchains = prepare_newgen(tmcmc_data, nchains, leaders);	// calculate statistics 

//	spmd_update_gdata();
//	spmd_print_matrix_2d();
	call_print_matrix_2d(tmcmc_data);

	/* this can be moved above */
	if (tmcmc_data->runinfo.p[tmcmc_data->runinfo.Gen] == 1) {
		printf("p == 1 from previous run, nothing more to do\n");
		goto end;
	}

#if 0
	for (runinfo.Gen = 1; tmcmc_data->runinfo.Gen < tmcmc_data->data.MaxStages; tmcmc_data->runinfo.Gen++){
#else
	tmcmc_data->runinfo.Gen++;
	for (/*runinfo.Gen = 1*/; tmcmc_data->runinfo.Gen < tmcmc_data->data.MaxStages; tmcmc_data->runinfo.Gen++){
#endif
		/* process current generation, compute probs, find new chains */
		//leader[i]: { point[tmcmc_data->data.Nth], F, nsteps}

		int dim = tmcmc_data->data.Nth;
		double bSS[dim*dim];

		for (i = 0; i < dim; i++) {
			int j;
			for (j = 0; j < dim; j++) {
				bSS[i*dim+j]= tmcmc_data->data.bbeta*tmcmc_data->runinfo.SS[i][j];
				printf("bSS[%d] = %f\n", i*dim+j, bSS[i*dim+j]/tmcmc_data->data.bbeta);
			}
		}



		int winfo[4];
		double in_tparam[tmcmc_data->data.Nth];
		int nsteps;
		gt0 = torc_gettime();
		
		for (i = 0; i < nchains; i++) {
			winfo[0] = tmcmc_data->runinfo.Gen;
			winfo[1] = i;
			winfo[2] = -1;	// not used
			winfo[3] = -1;	// not used

			int p;
			for (p = 0; p < tmcmc_data->data.Nth; p++)
				in_tparam[p] = leaders[i].point[p];
			nsteps = leaders[i].nsel;

			out_tparam[i] = leaders[i].F;	// fleader...

#ifdef USE_TORC
			torc_create(leaders[i].queue, chaintask, 11,
//			torc_create(torc_worker_id(), chaintask, 11,
				tmcmc_data->data.Nth, MPI_DOUBLE, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_REF,
				4, MPI_INT, CALL_BY_COP,
				1, MPI_DOUBLE, CALL_BY_COP,
				dim*dim, MPI_DOUBLE, CALL_BY_VAL,
				dim, MPI_DOUBLE, CALL_BY_VAL,
				dim, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_LONG, CALL_BY_VAD,
				2, MPI_INT, CALL_BY_COP,
				in_tparam, &tmcmc_data->data.Nth, &nsteps, &out_tparam[i], winfo,
				&tmcmc_data->runinfo.p[tmcmc_data->runinfo.Gen],		// pj
				bSS, tmcmc_data->data.lowerbound, tmcmc_data->data.upperbound,
				tmcmc_data, tmcmc_owner);
#else
			chaintask(
				in_tparam, &tmcmc_data->data.Nth, &nsteps, &out_tparam[i], winfo,
				&tmcmc_data->runinfo.p[tmcmc_data->runinfo.Gen],		// pj
				bSS, tmcmc_data->data.lowerbound, tmcmc_data->data.upperbound,
				tmcmc_data, tmcmc_owner);
#endif

		}
#ifdef USE_TORC
		/* wait for all chain tasks to finish */
#ifdef _STEALING_
		torc_enable_stealing();
#endif
		torc_waitall();
#ifdef _STEALING_
		torc_disable_stealing();
#endif
#endif

		gt1 = torc_gettime();
		int g_nfeval = get_nfc();
		printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", tmcmc_data->runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
		reset_nfc();

		if (tmcmc_data->data.idump) dump_curgen_db(&tmcmc_data->curgen_db, tmcmc_data->runinfo.Gen, tmcmc_data->data.Nth, tmcmc_info);





		// save here
		nchains = prepare_newgen(tmcmc_data, nchains, leaders);	// calculate statistics

		spmd_update_gdata();
		//spmd_print_matrix_2d();
		call_print_matrix_2d(tmcmc_data);

#if 0
		printf("=================\n");
		print_matrix("runinfo.p", tmcmc_data->runinfo.p, tmcmc_data->runinfo.Gen+1);
		print_matrix("runinfo.CoefVar", tmcmc_data->runinfo.CoefVar, tmcmc_data->runinfo.Gen+1);
		print_matrix_i("runinfo.currentuniques", tmcmc_data->runinfo.currentuniques, tmcmc_data->runinfo.Gen+1);
		print_matrix("runinfo.acceptance", tmcmc_data->runinfo.acceptance, tmcmc_data->runinfo.Gen+1);
		print_matrix("runinfo.logselection", tmcmc_data->runinfo.logselection, tmcmc_data->runinfo.Gen+1);
		printf("=================\n");
#endif

		if (tmcmc_data->runinfo.p[tmcmc_data->runinfo.Gen] == 1) {
			break;
		}
		if (tmcmc_data->runinfo.Gen+1 == tmcmc_data->data.MaxStages) {
			break;
		}
	}

	print_matrix("runinfo.p", tmcmc_data->runinfo.p, tmcmc_data->runinfo.Gen+1);
	print_matrix("runinfo.CoefVar", tmcmc_data->runinfo.CoefVar, tmcmc_data->runinfo.Gen+1);
	print_matrix_i("runinfo.currentuniques", tmcmc_data->runinfo.currentuniques, tmcmc_data->runinfo.Gen+1);
	print_matrix("runinfo.acceptance", tmcmc_data->runinfo.acceptance, tmcmc_data->runinfo.Gen+1);
	print_matrix("runinfo.logselection", tmcmc_data->runinfo.logselection, tmcmc_data->runinfo.Gen+1);

	double logEvidence[1];
	logEvidence[0] = compute_sum(tmcmc_data->runinfo.logselection, tmcmc_data->runinfo.Gen+1);
	print_matrix("logEvidence", logEvidence, 1);

	print_matrix_2d("runinfo.SS", tmcmc_data->runinfo.SS, tmcmc_data->data.Nth, tmcmc_data->data.Nth);

	for (i = 0; i < tmcmc_data->runinfo.Gen+1; i++) {
		char title[64];
		sprintf(title, "runinfo.meantheta(%d)", i);
		//print_matrix("runinfo.meantheta", tmcmc_data->runinfo.meantheta[i], tmcmc_data->data.Nth);
		print_matrix(title, tmcmc_data->runinfo.meantheta[i], tmcmc_data->data.Nth);
	}

	// last save here - do we need this? what happens if I restart the program with this saved data
end:
	/* terminate spawners */

	printf("total function calls = %d\n", get_tfc());

//	return logEvidence[0];
	*res = logEvidence[0];
}


#if 0
typedef void (*fitfun_t)();
fitfun_t fitfun_initialize;
fitfun_t fitfun_finalize;
#else
extern void fitfun_initialize(char *fitfun_name);
extern void fitfun_finalize();
extern double fitfun(double *TP, long n, void *output, int *info);
#endif

int tmcmc_initialize(char *fitfun_name)
{
	static int initialized = 0;

	if (initialized) return 0;
	initialized=1; 

#ifdef USE_TORC
	torc_register_task(initchaintask);
	torc_register_task(chaintask);
	torc_register_task(torc_update_curgen_db_task);
	torc_register_task(reset_nfc_task);
	torc_register_task(get_nfc_task);
	torc_register_task(taskfun);
	torc_register_task(call_gsl_rand_init);
	torc_register_task(call_print_matrix_2d);
	torc_register_task(call_update_gdata);
	torc_register_task(tmcmc);
#endif


#if 0
	void* handle = dlopen("libfitfun.so", RTLD_LAZY);
	if (handle == NULL)
	{
		printf("Error: Could not open libfitfun.so!\n");
		exit(1);
	}

	// reset errors
	dlerror();
	fitfun_initialize = (fitfun_t) dlsym(handle, "fitfun_initialize");
	fitfun_finalize = (fitfun_t) dlsym(handle, "fitfun_finalize");
	fitfun = (feval_t) dlsym(handle, "fitfun");
#endif


#ifdef USE_TORC
	torc_register_task(fitfun_initialize);
	torc_register_task(fitfun_finalize);

	torc_init(0, NULL, MODE_MS);

	int i;
	for (i = 0; i < torc_num_nodes(); i++)
		torc_create(i*torc_i_num_workers(), fitfun_initialize, 0);
	torc_waitall();
#else
	printf("calling fitfun_initialize (%s)\n", fitfun_name);
	fitfun_initialize(fitfun_name);
#endif
	return 0;
}

void tmcmc_finalize()
{
	static int finalized = 0;

	if (finalized) return;
	finalized=1; 

#ifdef USE_TORC
	int i;
	for (i = 0; i < torc_num_nodes(); i++)
		torc_create(i*torc_i_num_workers(), fitfun_finalize, 0);
	torc_waitall();

	torc_finalize();
#else
	fitfun_finalize();
#endif
}

