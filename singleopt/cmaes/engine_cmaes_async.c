/* --------------------------------------------------------- */
/* --------------- A Very Short Example -------------------- */
/* --------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h> /* free() */
#include "cmaes_interface.h"
#include <mpi.h>
#include <torc.h>
#include <unistd.h>
#include <math.h>

#define VERBOSE 1
//#define _STEALING_
//#define _RESTART
#define IODUMP 1
#include "fitfun.c" 


#if 1

#define MAXDIM 16
typedef struct data_s
{
	double in[MAXDIM];
	double out;
	int node_id;
	int state;
	int gen;	// generation
} data_t;

data_t *result;
int DATA_ENTRIES;

pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
int completed = 0;
int exit_flag = 0;

int find_tid()
{
	int i;
	for (i = 0; i < DATA_ENTRIES; i++)
	{
		if (result[i].state == 0)
		{
			result[i].state = 1;
			return i;
		}
	}

	printf("This should not happen!\n");
	exit(1);
	return -1;
}

int current_jobs()
{
	int i;

	int s = 0;
	for (i = 0; i < DATA_ENTRIES; i++)
	{
		if ((result[i].state == 1) || (result[i].state == 2))
		{
			s++;
		}
	}
	return s;
}


void setexit()
{
	exit_flag = 1;
}

int extra;
int lambda, dim;
int step = 0;
int taskid = 0; 

cmaes_t evo; /* an CMA-ES type struct or "object" */
double *arFunvals, *const*pop, *xfinal;

int is_feasible(double *pop, int dim);
void taskfun(int *ptid, double *x, int *pn, int *info);


void callback(int *ptid, int *pnode_id, double *pout)
{
	int tid = *ptid;
	int node_id = *pnode_id;
	double out = *pout;

	pthread_mutex_lock(&m);
	completed++;

	printf("vvvvvvvvvvvvvvvv\n"); fflush(0);
	printf("cb: tid = %d, node_id= %d pout = %f, completed = %d\n", tid, node_id, out, completed); fflush(0);
	result[tid].out = out;
	result[tid].node_id = node_id;
	result[tid].state = 2;  // done

	if (completed % lambda == 0)
	{
		printf("!!! completed == %d\n", completed); fflush(0);
		printf("^^^^^^^^^^^^^^^^\n"); fflush(0);

		int i;
		int j = 0;
		for (i = 0; i < DATA_ENTRIES; i++)
		{
			if (result[i].state == 2) {
				int k;
				for (k = 0; k < dim; k++) pop[j][k] = result[i].in[k];
				evo.publicFitness[j] = result[i].out;
				arFunvals[j] = result[i].out;	//pop[j][dim] = result[i].out;
				//printf("setting arFunvals[%d]=%f to result[%d].out=%f\n", j, arFunvals[j], i, result[i].out); 
				j++;
			}
		}

#if 0
		printf("arFunvals = [\n");
		for (i = 0; i < lambda; i++)
		{
			printf("%f\n", arFunvals[i]);
		}
		printf("];\n");
#endif

		cmaes_UpdateDistribution(&evo, arFunvals);  

		/* generate lambda new search points, sample population */
		/* read instructions for printing output or changing termination conditions */ 
		cmaes_ReadSignals(&evo, "cmaes_signals.par");   
		fflush(stdout); /* useful in MinGW */

#if VERBOSE
		{
		const double *xbever = cmaes_GetPtr(&evo, "xbestever");
		double fbever = cmaes_Get(&evo, "fbestever");
	
		printf("BEST @ %5d: ", step);
		for (i = 0; i < dim; i++)
			printf("%25.16lf ", xbever[i]);
		printf("%25.16lf\n", fbever);
		}
//		printf("Step %4d: time = %.3lf seconds\n", step, tt1-tt0);
#endif

#if IODUMP
		{
		char filename[256];
		sprintf(filename, "curgen_db_%03d.txt", step);
		FILE *fp = fopen(filename, "w");
		for (i = 0; i < lambda; i++) {
			int j;
			for (j = 0; j < dim; j++) fprintf(fp, "%.6e ", pop[i][j]);
			fprintf(fp, "%.6e\n", arFunvals[i]);
		}
		fclose(fp);
		}
#endif


#if defined(_RESTART_)
		cmaes_WriteToFile(&evo, "resume", "allresumes.dat");         /* write resume data */
#endif


		if (cmaes_TestForTermination(&evo))
		{
			exit_flag = 1;
			pthread_mutex_unlock(&m);
			return;
		}

		for (i = 0; i < DATA_ENTRIES; i++)
		{
			//if (result[i].state == 2) {
				printf("%d: %f %f %f %d\n", i, result[i].in[0], result[i].in[1], result[i].out, result[i].state);
			//}
			if (result[i].state == 2) {
				result[i].state = 0;    // free
				result[i].in[0]=-1;
				result[i].in[1]=-1;
				result[i].out=-1;
			}
		}

		step++;
		taskid = 0;	// reset

		pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */
		for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i)
			while (!is_feasible(pop[i], dim))
				cmaes_ReSampleSingle(&evo, i);


		printf("tasksize = %d+%d = %d\n", lambda, extra, lambda + extra);
		int tospawn = lambda;

#if 0
		if ((step >= 40) && (step % 10 == 0) && (extra > 0))
		{
			extra--;
			tospawn--;
		}
#endif


		for (i = 0; i < tospawn; i++)
		{
			int tid = find_tid();
			int j;
			for (j = 0; j < dim; j++) result[tid].in[j] = pop[i][j];

	                printf("yy: spawning task %d\n", tid);
			
			int info[4];
			info[0] = 0; info[1] = 0; info[2] = step; info[3] = taskid++; 	/* gen, chain, step, task */
			torc_create_detached(tid % torc_num_workers(), taskfun, 4,
				1, MPI_INT, CALL_BY_COP,
				dim, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_INT, CALL_BY_COP,
				4, MPI_INT, CALL_BY_COP,
				&tid, result[tid].in, &dim, info);
		}
	}
	else
	{
		if (exit_flag) 
		{
			pthread_mutex_unlock(&m);
			return;
		}

		int njobs = current_jobs();
		if (njobs < (lambda + extra))
		{
			int tid = find_tid();
			int i = tid % lambda;	// could be zero

			cmaes_ReSampleSingle(&evo, i);
			while (!is_feasible(pop[i], dim)) 
			{
				cmaes_ReSampleSingle(&evo, i); 
			}

			int j;
			for (j = 0; j < dim; j++) result[tid].in[j] = pop[i][j];

			printf("yy: spawning task %d\n", tid);
			
			int info[4];
			info[0] = 0; info[1] = 0; info[2] = step; info[3] = taskid++; 	/* gen, chain, step, task */
			torc_create_detached(tid % torc_num_workers(), taskfun, 4,
				1, MPI_INT, CALL_BY_COP,
				dim, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_INT, CALL_BY_COP,
				4, MPI_INT, CALL_BY_COP,
				&tid, result[tid].in, &dim, info);
		}

	}

	pthread_mutex_unlock(&m);
	return;
}

extern void _torc_scheduler_loop(int);
void torc_dispatch()
{
        _torc_scheduler_loop(1);
}

#endif


double *lower_bound;	//double lower_bound[] = {-6.0, -6.0};
double *upper_bound;	//double upper_bound[] = {+6.0, +6.0};


/* the objective (fitness) function to be minimized */
void taskfun(int *ptid, double *x, int *pn, int *info)
{
	int tid = *ptid;
	int n = *pn;
	int gen, chain, step, task;
	gen = info[0]; chain = info[1]; step = info[2]; task = info[3];
	printf("executing task (%d,%d,%d,%d)\n", gen, chain, step, task);
	
	double f = -fitfun(x, n, (void *)NULL, info);	/* CMA-ES needs this minus sign */

	printf("taskfun %d: %f %f = %f\n", tid, x[0], x[1], f);

	int node_id = torc_node_id();

	if (node_id == 0)
	{
		callback(&tid, &node_id, &f);
	}
	else
	{
		torc_create_direct(0, callback, 3,
			1, MPI_INT, CALL_BY_COP,
			1, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_COP,
			&tid, &node_id, &f);
		torc_waitall3();
	}

	return;
}


int is_feasible(double *pop, int dim)
{
	int i, good;
	// printf("is_feasible %d\n", dim);
	for (i = 0; i < dim; i++) {
		good = (lower_bound[i] <= pop[i]) && (pop[i] <= upper_bound[i]);
		if (!good) {
			//printf("%d not good\n", i);
			//usleep(1000);
			return 0;
		}
	}
	return 1;
}


/* the optimization loop */
int main(int argn, char **args)
{
	int i; 
	/* peh - start */
	//int lambda, dim;  - peh: global variables
	//int step = 0; - peh: global variable
	double gt0, gt1, gt2, gt3;
	//double tt0 = 0, tt1 = 0, stt = 0.0; 
	int info[4];	/* gen, chain, step, task */
	/* peh - end */

	torc_register_task(taskfun);
	torc_register_task(callback);

	/* Initialize everything into the struct evo, 0 means default */
	torc_init(argn, args, MODE_MS);

	gt0 = torc_gettime();
	arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "cmaes_initials.par"); 

	printf("%s\n", cmaes_SayHello(&evo));

	cmaes_ReadSignals(&evo, "cmaes_signals.par");  /* write header and initial values */

	lambda = cmaes_Get(&evo, "lambda"); 
	dim = cmaes_Get(&evo, "dim"); 

	lower_bound = malloc(dim*sizeof(double));
	upper_bound = malloc(dim*sizeof(double));
	for (i = 0; i < dim; i++) {
		lower_bound[i] = -6;
		upper_bound[i] = +6;
	}

	gt1 = torc_gettime();

	/* generate lambda new search points, sample population */
	pop = cmaes_SamplePopulation(&evo); /* do not change content of pop */
	for (i = 0; i < cmaes_Get(&evo, "popsize"); ++i) 
		while (!is_feasible(pop[i], dim)) 
			cmaes_ReSampleSingle(&evo, i); 


	extra = lambda - 1;
	int cnt = lambda + extra;
	DATA_ENTRIES = cnt; // + lambda;

	result = (data_t *)calloc(1, DATA_ENTRIES*sizeof(data_t));
	for (i = 0; i < lambda; i++)
	{
		int j;
		for (j = 0; j < dim; j++) result[i].in[j] = pop[i][j];
		result[i].state = 0;
	}


	for (i=0; i<lambda; i++)
	{
		pthread_mutex_lock(&m);
		int tid = find_tid();
		pthread_mutex_unlock(&m);

                printf("yy: spawning task %d\n", tid);

		info[0] = 0; info[1] = 0; info[2] = step; info[3] = taskid++; 	/* gen, chain, step, task */
		torc_create_detached(tid % torc_num_workers(), taskfun, 4,
				1, MPI_INT, CALL_BY_COP,
				dim, MPI_DOUBLE, CALL_BY_VAL,
				1, MPI_INT, CALL_BY_COP,
				4, MPI_INT, CALL_BY_COP,
				&tid, result[tid].in, &dim, info);
        }

        while (!exit_flag)
        {
                torc_dispatch();
        }

	gt2 = torc_gettime();


	printf("Stop:\n%s\n",  cmaes_TestForTermination(&evo)); /* print termination reason */
	cmaes_WriteToFile(&evo, "all", "allcmaes.dat");         /* write final results */

	/* get best estimator for the optimum, xmean */
	xfinal = cmaes_GetNew(&evo, "xmean"); /* "xbestever" might be used as well */
	cmaes_exit(&evo); /* release memory */ 

	/* do something with final solution and finally release memory */
	free(xfinal); 

	gt3 = torc_gettime();
	printf("Total elapsed time = %.3lf  seconds\n", gt3-gt0);
	printf("Initialization time = %.3lf  seconds\n", gt1-gt0);
	printf("Processing time = %.3lf  seconds\n", gt2-gt1);
	//printf("Funtion Evaluation time = %.3lf  seconds\n", stt);
	printf("Finalization time = %.3lf  seconds\n", gt3-gt2);

	torc_finalize();

	return 0;
}

