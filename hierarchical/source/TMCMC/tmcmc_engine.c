/*
 *  engine_tmcmc.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <assert.h>
#include <signal.h>
#include <stdio.h>
#include <unistd.h>

#include "tmcmc_db.h"
#include "tmcmc_aux.h"
#include "tmcmc_stats.h"
#include "tmcmc_engine.h"
#include "fitfun.h"


#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


#include "../priors/priors.h"
#include "../priors/myrand.h"
Density *priors;


// #define _STEALING_
// #define VERBOSE 1
// #define _RESTART_






data_t 		data;
runinfo_t 	runinfo;
cgdb_t 		curgen_db;
db_t 		full_db;
resdb_t 	curres_db;







void data_init()
{
    int i;

    /* DATA: user's input parameters */
    read_data();

    init_curgen_db();
    init_curres_db();
    init_full_db();

    /* RUNINFO: running state */
    runinfo.CoefVar 		= (double *)calloc(1, (data.MaxStages+1)*sizeof(double));
    runinfo.p 				= (double *)calloc(1, (data.MaxStages+1)*sizeof(double));
    runinfo.currentuniques 	= (int 	  *)calloc(1, data.MaxStages*sizeof(int));
    runinfo.logselection 	= (double *)calloc(1, data.MaxStages*sizeof(double));
    runinfo.acceptance 		= (double *)calloc(1, data.MaxStages*sizeof(double));

    double *SSmem = (double *)calloc(1, data.Nth*data.Nth*sizeof(double));
    
	runinfo.SS    = (double **)malloc(data.Nth*sizeof(double *));
	for (i = 0; i < data.Nth; i++){
        runinfo.SS[i] = SSmem + i*data.Nth; /*&SSmem[i*data.Nth];*/
    }

    runinfo.meantheta = (double **)calloc(1, (data.MaxStages+1)*sizeof(double *));
    for (i = 0; i < data.MaxStages+1; i++) {
        runinfo.meantheta[i] = (double *)calloc(1, data.Nth*sizeof(double));
    }

    runinfo.Gen = 0;
    runinfo.CoefVar[0] = 10;

	if(data.options.Display){
    	printf("runinfo    = %p\n", &runinfo);
    	printf("runinfo.p  = %p\n", runinfo.p);
    	printf("runinfo.SS = %p\n", runinfo.SS);
	}

}






void read_data()
{
    int i;

    data.Nth 		= -1;
    data.MaxStages 	= -1;
    data.PopSize 	= -1;

    data.MinChainLength = 0;
    data.MaxChainLength = 1; //1e6;

    data.TolCOV 	= -1;
    data.bbeta  	= -1;
    data.seed 		= 280675;
    data.burn_in 	= -1;

    data.options.MaxIter 	= -1; 
    data.options.Tol 		= -1;
    data.options.Display 	= 0;
    data.options.Step 		= -1;

    data.load_from_file = 0;    /* load initial samples from file instead from prior */

    data.icdump = 1;    /* dump current dataset of accepted points */
    data.ifdump = 0;    /* dump complete dataset of points */

    data.stealing = 0;
    data.restart = 0;



    // Read parameters from file
    FILE *f = fopen("tmcmc.par", "r");
    if (f == NULL) {
		printf("\nThe input file 'tmcmc.par' is missing. Exit...\n");
		exit(1);
    }

    char line[256];

    int line_no = 0;
    while (fgets(line, 256, f)!= NULL) {
        line_no++;
        if ((line[0] == '#')||(strlen(line)==0)) {
            continue;
        }

        if (strstr(line, "Nth")) {
            sscanf(line, "%*s %d", &data.Nth);
        }
        else if (strstr(line, "MaxStages")) {
            sscanf(line, "%*s %d", &data.MaxStages);
        }
        else if (strstr(line, "PopSize")) {
            sscanf(line, "%*s %d", &data.PopSize);
        }
        else if (strstr(line, "TolCOV")) {
            sscanf(line, "%*s %lf", &data.TolCOV);
        }
        else if (strstr(line, "bbeta")) {
            sscanf(line, "%*s %lf", &data.bbeta);
        }
        else if (strstr(line, "seed")) {
            sscanf(line, "%*s %ld", &data.seed);
        }
        else if (strstr(line, "burn_in")) {
            sscanf(line, "%*s %d", &data.burn_in);
        }
        else if (strstr(line, "opt.MaxIter")) {
            sscanf(line, "%*s %d", &data.options.MaxIter);
        }
        else if (strstr(line, "opt.Tol")) {
            sscanf(line, "%*s %lf", &data.options.Tol);
        }
        else if (strstr(line, "opt.Display")) {
            sscanf(line, "%*s %d", &data.options.Display);
        }
        else if (strstr(line, "opt.Step")) {
            sscanf(line, "%*s %lf", &data.options.Step);
            printf("setting step = %f\n", data.options.Step);
        }
        else if (strstr(line, "icdump")) {
            sscanf(line, "%*s %d", &data.icdump);
        }
        else if (strstr(line, "ifdump")) {
            sscanf(line, "%*s %d", &data.ifdump);
        }
        else if (strstr(line, "use_local_cov")) {
            sscanf(line, "%*s %d", &data.use_local_cov);
        }
        else if (strstr(line, "stealing")) {
            sscanf(line, "%*s %d", &data.stealing);
        }
        else if (strstr(line, "restart")) {
            sscanf(line, "%*s %d", &data.restart);
        }
    }

	fclose(f);


	//XXX add: check if all parameters are ok
	

    data.lowerbound = (double *)malloc(data.Nth*sizeof(double));
    data.upperbound = (double *)malloc(data.Nth*sizeof(double));
	for(i=0; i<data.Nth; i++){
 		data.lowerbound[i] = -INFINITY;
    	data.upperbound[i] =  INFINITY;
	}



    data.Num = (int *)malloc(data.MaxStages*sizeof(int));
    for (i = 0; i < data.MaxStages; i++){
        data.Num[i] = data.PopSize;
    }
    data.LastNum = data.PopSize;

    
	
	double *LCmem  = (double  *)calloc(1, data.PopSize*data.Nth*data.Nth*sizeof(double));
    data.local_cov = (double **)malloc(data.PopSize*sizeof(double *));
    for (int pos=0; pos < data.PopSize; ++pos){
        data.local_cov[pos] = LCmem + pos*data.Nth*data.Nth;
        for (i=0; i<data.Nth; ++i)
            data.local_cov[pos][i*data.Nth+i] = 1;
    }


}







void save_runinfo()
{
    int i, j;

    /* allocate and initialize runinfo */
    FILE *f = fopen("runinfo.txt", "w");

    fprintf(f, "Gen=\n");
    fprintf(f, "%d\n", runinfo.Gen);

    fprintf(f, "CoefVar[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.CoefVar[i]);

    fprintf(f, "p[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.p[i]);

    fprintf(f, "currentuniques[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%d\n", runinfo.currentuniques[i]);

    fprintf(f, "logselection[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.logselection[i]);

    fprintf(f, "acceptance[%d]=\n", data.MaxStages);
    for (i = 0; i < data.MaxStages; i++) fprintf(f, "%.16lf\n", runinfo.acceptance[i]);

    fprintf(f, "SS[%d][%d]=\n", data.Nth, data.Nth);
    for (i = 0; i < data.Nth; i++)
        for (j = 0; j < data.Nth; j++)
            fprintf(f, "%.16lf\n", runinfo.SS[i][j]);

    fprintf(f, "meantheta[%d][%d]\n", data.MaxStages, data.Nth);
    for (i = 0; i < data.MaxStages; i++)
        for (j = 0; j < data.Nth; j++)
            fprintf(f, "%.16lf\n", runinfo.meantheta[i][j]);

    fclose(f);
}







int load_runinfo()
{
    int i, j;
    char header[256];

    /* allocate and initialize runinfo */
    FILE *f = fopen("runinfo.txt", "r");
    if (f == NULL) return 1;

    fscanf(f, "%s", header);
    fscanf(f, "%d", &runinfo.Gen);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.CoefVar[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.p[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%d\n", &runinfo.currentuniques[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.logselection[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++) fscanf(f, "%lf\n", &runinfo.acceptance[i]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.Nth; i++)
        for (j = 0; j < data.Nth; j++)
            fscanf(f, "%lf\n", &runinfo.SS[i][j]);

    fscanf(f, "%s", header);
    for (i = 0; i < data.MaxStages; i++)
        for (j = 0; j < data.Nth; j++)
            fscanf(f, "%lf\n", &runinfo.meantheta[i][j]);

    fclose(f);

    return 0;
}










void initchaintask(double in_tparam[], int *pdim, double *out_tparam, int winfo[4])
{
    int i;
    int gen_id, chain_id;
    gen_id = winfo[0];
    chain_id = winfo[1];

    long 	me = torc_worker_id();
    double 	point[data.Nth], fpoint;


	//printf("-->%ld / %d\n", me, torc_i_num_workers() );

    for (i = 0; i < data.Nth; i++)
        point[i] = in_tparam[i];

    evaluate_F(point, &fpoint, me, gen_id, chain_id, 0, 1);


	double logprior = prior_log_pdf( priors, data.Nth, point );

    /* update current db entry */
    torc_update_curgen_db( point, fpoint, logprior );
    if (data.ifdump) torc_update_full_db(point, fpoint, NULL, 0, 0);
    *out_tparam = fpoint;    /* currently not required, the result is already in the db*/


    return;
}




void evaluate_F(double point[], double *Fval, int worker_id, int gen_id, int chain_id, int step_id, int ntasks)
{
    double F;
    int winfo[4];
    int dim = data.Nth;

    winfo[0] = gen_id;
    winfo[1] = chain_id;
    winfo[2] = step_id;
    winfo[3] = 0;

	#if VERBOSE
    	printf("running on worker %d\n", worker_id);
	#endif
    
	taskfun(point, &dim, &F, winfo);

    *Fval = F;
}





void taskfun(double /*const*/ *x, int *pN, double *res, int winfo[4])
{
    double f;
    int N = *pN;

    inc_nfc();    /* increment function call counter*/

    f = fitfun(x, N, (void *)NULL, winfo);
	#if (EXPERIMENTAL_RESULTS > 0)    /* peh: decide about this (results should be passed as argument to fitfun) */
		double results[EXPERIMENTAL_RESULTS];
		for (int i = 0; i < EXPERIMENTAL_RESULTS; i++) {
			if (i < data.Nth)
				results[i] = x[i];
			else
				results[i] = 0.0;
		}
		torc_update_curres_db(results, f);
	#endif

    *res = f;
    return;
}







//------------------------------------------------------------------------
//		Update database functions
//

void torc_update_curgen_db(double point[], double F, double prior)
{
    int me = torc_node_id();

    if (me == 0) {
        update_curgen_db(point, F,prior);
        return;
    }

	#if defined(_USE_TORC_)
		torc_create_direct(0, (void (*)())torc_update_curgen_db_task, 3,
						/* message to the database manager (separate process?) or direct execution by server thread */
			data.Nth, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_COP,
			point, &F,&prior);

    		torc_waitall3();    /* wait without releasing the worker */
	#endif

}





void torc_update_full_db(double point[], double F, double *G, int n, int surrogate)
{
    if (torc_node_id() ==0) {
        update_full_db(point, F, G, n, surrogate);
        return;
    }

	#if defined(_USE_TORC_)
    	if (n == 0)
        	torc_create_direct( 0, (void (*)())torc_update_full_db_task, 3,        
									/* message to the database manager (separate process?) or direct execution by server thread */
                				data.Nth, MPI_DOUBLE, CALL_BY_VAL,
                				1, MPI_DOUBLE, CALL_BY_COP,
                				1, MPI_INT, CALL_BY_COP,
                				point, &F, &surrogate);
    	else
        	torc_create_direct( 0, (void (*)())torc_update_full_db_task_p5, 5,
                					data.Nth, MPI_DOUBLE, CALL_BY_VAL,
                					1, MPI_DOUBLE, CALL_BY_COP,
                					n, MPI_DOUBLE, CALL_BY_COP,    /* xxx: for CALL_BY_VAL: in the full-version of the library, with n=1 we had segv */
                					1, MPI_INT, CALL_BY_COP,
                					1, MPI_INT, CALL_BY_COP,
                					point, &F, G, &n, &surrogate);

    		torc_waitall3();
	#endif
}







void torc_update_full_db_task(double point[], double *pF, int *psurrogate)
{
    double F = *pF;
    int surrogate = *psurrogate;
    double *G = NULL;
    int n = 0;

    update_full_db(point, F, G, n, surrogate);
}



void torc_update_full_db_task_p5(double point[], double *pF, double *G, int *pn, int *psurrogate)
{
    double F = *pF;
    int n = *pn;
    int surrogate = *psurrogate;

    update_full_db(point, F, G, n, surrogate);
}






void torc_update_curgen_db_task(double point[], double *pF, double *pprior)
{
    double F = *pF;
	double prior = *pprior;

    update_curgen_db(point, F, prior);
}














//------------------------------------------------------------------------
//		Exit functions
//
void check_for_exit()
{
    int val, exitgen = -1;
    char *s;

    s = (char *) getenv("EXIT_GEN");
    if (s != 0 && sscanf(s, "%d", &val) == 1 && val >= 0)
        exitgen = val;

    if (exitgen == runinfo.Gen) {
        printf("Read Exit Envrironment Variable!!!\n");

		#if defined(_USE_TORC_)
        	torc_finalize();
		#endif
        exit(1);
    }

    FILE *fp;
    fp = fopen("exit.txt", "r");
    if (fp != NULL) {
        printf("Found Exit File!!!\n");
        unlink("exit.txt");
		#if defined(_USE_TORC_)
        	torc_finalize();
		#endif
        exit(1);
    }
}






















void torc_update_curres_db_task(double point[EXPERIMENTAL_RESULTS], double *pF)
{
    double F = *pF;

    update_curres_db(point, F);
}








void torc_update_curres_db(double point[EXPERIMENTAL_RESULTS], double F)
{
    int me = torc_node_id();

    if (me ==0) {
        update_curres_db(point, F);
        return;
    }

	#if defined(_USE_TORC_)
    	torc_create_direct(	0, (void (*)())torc_update_curres_db_task, 2,        
							/* message to the database manager (separate process?) or direct execution by server thread */
            				EXPERIMENTAL_RESULTS, MPI_DOUBLE, CALL_BY_COP,
            				1, MPI_DOUBLE, CALL_BY_COP,
            				point, &F);
    	torc_waitall3();
	#endif
}






double F(double *TP, int *pn)    /* for PNDL */
{
    double gres;

    taskfun(TP, pn, &gres, NULL);

    return gres;
}












static int in_rect(double *v1, double *v2, double *diam, double sc, int D)
{
    int d;
    for (d = 0; d < D; ++d) {
        if (fabs(v1[d]-v2[d]) > sc*diam[d]) return 0;
    }
    return 1;
}






void precompute_chain_covariances(const cgdbp_t* leader,double** init_mean, double** chain_cov, int newchains)
{
    printf("Precomputing covariances for the current generation...\n");

    int i, j, k, d, pos, ind;

    int D = data.Nth;
    int N = curgen_db.entries;

    double my_time = clock();

    // allocate space
    int* nn_ind = (int*)malloc(newchains*N*sizeof(int));
    int* nn_count = (int*)malloc(newchains*sizeof(int));
    double* diam = (double*)malloc(D*sizeof(double));
    double* chain_mean = (double*)malloc(D*sizeof(double));
    gsl_matrix* work = gsl_matrix_alloc(D, D);

    // find diameters
    for (d = 0; d < D; ++d) {
        double d_min = +1e6;
        double d_max = -1e6;
        for (pos = 0; pos < N; ++pos) {
            double s = curgen_db.entry[pos].point[d];
            if (d_min > s) d_min = s;
            if (d_max < s) d_max = s;
        }
        diam[d] = d_max-d_min;
        printf("Diameter %d: %.6lf\n", d, diam[d]);
    }

    int status = 0;
    double scale, ds = 0.05;
    for (scale = 0.1; scale <= 1.0; scale += ds) {
        // find neighbors in a rectangle - O(N^2)
        for (pos = 0; pos < newchains; ++pos) {
            nn_count[pos] = 0;
            double* curr = leader[pos].point;
            for (i = 0; i < N; ++i) {
                double* s = curgen_db.entry[i].point;
                if (in_rect(curr, s, diam, scale, D)) {
                    nn_ind[pos*N+nn_count[pos]] = i;
                    nn_count[pos]++;
                }
            }
        }

        // compute the covariances
        for (pos = 0; pos < newchains; ++pos) {
            for (d = 0; d < D; ++d) {
                chain_mean[d] = 0;
                for (k = 0; k < nn_count[pos]; ++k) {
                    ind = nn_ind[pos*N+k];
                    chain_mean[d] += curgen_db.entry[ind].point[d];
                }
                chain_mean[d] /= nn_count[pos];
            }

            for (i = 0; i < D; i++)
                for (j = 0; j < D; j++) {
                    double s = 0;
                    for (k = 0; k < nn_count[pos]; k++) {
                        ind = nn_ind[pos*N+k];
                        s += (curgen_db.entry[ind].point[i]-chain_mean[i]) *
                            (curgen_db.entry[ind].point[j]-chain_mean[j]);
                    }
                    chain_cov[pos][i*D+j] = chain_cov[pos][j*D+i] = s/nn_count[pos];
                }

            // check if the matrix is positive definite
            for (i = 0; i < D; ++i)
                for (j = 0; j < D; ++j){
                    double s = chain_cov[pos][i*D+j];
                    gsl_matrix_set(work, i, j, s);
                }
            gsl_set_error_handler_off();
            status = gsl_linalg_cholesky_decomp(work);
            if (status == GSL_SUCCESS) break;
        }
    }

    if (status != GSL_SUCCESS) {
        for (i = 0; i < D; i++)
            for (j = 0; j < D; j++)
                chain_cov[pos][i*D+j] = data.bbeta*runinfo.SS[i][j];
    }

    free(nn_ind);
    free(nn_count);
    free(diam);
    free(chain_mean);
    gsl_matrix_free(work);


    my_time = (clock() - my_time) / CLOCKS_PER_SEC;
    printf("Covariance computation time: %.2lf sec\n", my_time);
}





int compute_candidate(double candidate[], double chain_mean[], double var)
{
    int i, j;
    double bSS[data.Nth*data.Nth];

    for (i = 0; i < data.Nth; i++)
        for (j = 0; j < data.Nth; j++)
            bSS[i*data.Nth+j]= data.bbeta*runinfo.SS[i][j];

    
	mvnrnd(chain_mean, (double *)bSS, candidate, data.Nth);
    
	for (i = 0; i < data.Nth; i++) {
        if (isnan(candidate[i])) {
            printf("!!!!  isnan in candidate point!\n");
            exit(1);
            break;
        }
        if ((candidate[i] < data.lowerbound[i])||(candidate[i] > data.upperbound[i])) break;
    }
    
	if (i < data.Nth) return -1;

    return 0;    // all good
}







int compute_candidate_cov(double candidate[], double chain_mean[], double chain_cov[])
{
    int i;

    mvnrnd(chain_mean, (double *)chain_cov, candidate, data.Nth);
    for (i = 0; i < data.Nth; i++) {
        if (isnan(candidate[i])) {
            printf("!!!!  isnan in candidate point!\n");
            exit(1);
            break;
        }
        if ((candidate[i] < data.lowerbound[i])||(candidate[i] > data.upperbound[i])) return -1;
    }
    return 0;
}







void chaintask(double in_tparam[], int *pdim, int *pnsteps, double *out_tparam, int winfo[4],
        double *init_mean, double *chain_cov)
{

	

    int i,step;
    int nsteps = *pnsteps;
    int gen_id = winfo[0];
    int chain_id = winfo[1];

    long me = torc_worker_id();


    double leader[data.Nth], loglik_leader;        /* old*/
    double candidate[data.Nth], loglik_candidate;    /* new*/

    for (i = 0; i < data.Nth; i++)
        leader[i] = in_tparam[i];    /*chainwork->in_tparam[i];*/    /* get initial leader */
    loglik_leader = *out_tparam;    /*chainwork->out_tparam[0];*/    /* and its value */


    double pj = runinfo.p[runinfo.Gen];

    int burn_in = data.burn_in;

    for (step = 0; step < nsteps + burn_in; step++) {
        double chain_mean[data.Nth];
        if (step == 0)
            for (i = 0; i < data.Nth; ++i) chain_mean[i] = init_mean[i];
        else
            for (i = 0; i < data.Nth; ++i) chain_mean[i] = leader[i];

		//printf("---> %d -  %ld/%d   ---  %p  \n", chain_id, me, torc_i_num_workers(), priors ); fflush(NULL);

		#if 0
        	int fail = compute_candidate_cov(candidate, chain_mean, chain_cov);
		#else
        	int fail = compute_candidate(candidate, chain_mean, 1); // I keep this for the moment, for performance reasons
		#endif

        if (!fail){
            
			evaluate_F(candidate, &loglik_candidate, me, gen_id, chain_id, step, 1);    // this can spawn many tasks

            if (data.ifdump && step >= burn_in) torc_update_full_db(candidate, loglik_candidate, NULL, 0, 0);   
												// last argument should be 1 if it is a surrogate

            // Decide
            double logprior_candidate = prior_log_pdf( priors, data.Nth, candidate );
            double logprior_leader    = prior_log_pdf( priors, data.Nth, leader );
            
			double L;
            L = exp((logprior_candidate-logprior_leader)+(loglik_candidate-loglik_leader)*pj);

            if (L > 1) L = 1;
            double P = uniformrand(0,1);
            if (P < L) {
                for (i = 0; i < data.Nth; i++) leader[i] = candidate[i];    // new leader!
                loglik_leader = loglik_candidate;
                if (step >= burn_in) {
					double logprior_leader = prior_log_pdf( priors, data.Nth, leader );
					torc_update_curgen_db(leader, loglik_leader, logprior_leader);
				}
            }
            else {
                // increase counter or add the leader again in curgen_db
                if (step >= burn_in) {
					double logprior_leader = prior_log_pdf( priors, data.Nth, leader );
					torc_update_curgen_db(leader, loglik_leader, logprior_leader);
				}

            }
        }
		else{
            // increase counter or add the leader again in curgen_db
            if (step >= burn_in){
				double logprior_leader = prior_log_pdf( priors, data.Nth, leader );
				torc_update_curgen_db(leader, loglik_leader, logprior_leader);
			}
        }
    }

    return;
}











int compar_desc(const void* p1, const void* p2)
{
    int dir = +1;   /* -1: ascending order, +1: descending order */
    sort_t *s1 = (sort_t *) p1;
    sort_t *s2 = (sort_t *) p2;

    if (s1->nsel < s2->nsel) return dir;
    if (s1->nsel > s2->nsel) return -dir;
    /*    if (s1->nsel == s2->nsel) return 0;*/
    return 0;
}














int prepare_newgen(int nchains, cgdbp_t *leaders)
{
    /* process curgen_db -> calculate statitics */
    /* compute probs based on F values */
    /* draw new samples (nchains or user-specified) */
    /* find unique samples: fill the (new) leaders table */
    /* count how many times they appear -> nsteps */
    /* return the new sample size (number of chains) */

    int i, p;

    int n = curgen_db.entries;

    double *fj 		   = (double *) malloc( n*sizeof(double) );
    unsigned int *sel  = (unsigned int *) malloc( n*sizeof(sel) );

    double **g_x;
    g_x = (double **)malloc(data.Nth*sizeof(double *));
    for (i = 0; i < data.Nth; i++)
        g_x[i] = (double *)malloc(n*sizeof(double));

   
	if(0)	
	{
        double **x = g_x;

        for (p = 0; p < data.Nth; p++) {
            for (i = 0; i < n; i++) {
                x[p][i] = curgen_db.entry[i].point[p];
            }
        }

        double meanx[data.Nth], stdx[data.Nth];
        for (p = 0; p < data.Nth; p++) {
            meanx[p] = compute_mean(x[p], n);
            stdx[p] = compute_std(x[p], n, meanx[p]);
        }

		if(data.options.Display){
			printf("CURGEN DB (COMPLE) %d\n", runinfo.Gen);
			print_matrix((char *)"means", meanx, data.Nth);
	        print_matrix((char *)"std", stdx, data.Nth);
		}

    }



	if (1)
    {
        double **x = g_x;
        int un = 0, unflag, j;

        for( p = 0; p < data.Nth; p++ )
            x[p][un] = curgen_db.entry[0].point[p];
        
        un++;
        for (i = 1; i < n; i++){
            double xi[data.Nth];
            for (p = 0; p < data.Nth; p++)
            	xi[p] = curgen_db.entry[i].point[p];
           
            unflag = 1;    /* is this point unique?*/
            for (j = 0; j < un; j++){    /* check all the previous unique points*/
                int compflag;
                compflag = 1; 
                for (p = 0; p < data.Nth; p++){
                    if (fabs(xi[p]-x[p][j]) > 1e-6) {
                        compflag = 0;    /* they differ*/
                        break;
                    }
                }

                if (compflag == 1){
                    unflag = 0;    /* not unique, just found it in the unique points table*/
                    break;
                }
            }

            if (unflag){    /* unique, put it in the table */
                for (p = 0; p < data.Nth; p++) {
                    x[p][un] = xi[p];
                }
                un++;
            }
        }

        runinfo.currentuniques[runinfo.Gen] = un; /*+ 1*/;
        runinfo.acceptance[runinfo.Gen] = (1.0*runinfo.currentuniques[runinfo.Gen])/data.Num[runinfo.Gen]; /* check this*/

        double meanx[data.Nth], stdx[data.Nth];
        for (p = 0; p < data.Nth; p++) {
            meanx[p] = compute_mean(x[p], un);
            stdx[p] = compute_std(x[p], un, meanx[p]);
        }

        printf("CURGEN DB (UNIQUE) %d: [un = %d]\n", runinfo.Gen, un); /* + 1);*/
		if(data.options.Display){
			print_matrix((char *)"means", meanx, data.Nth);
	        print_matrix((char *)"std", stdx, data.Nth);
		}

    } /* end block*/






    {
		double t0 = torc_gettime();
		for (i = 0; i < n; i++)
			fj[i] = curgen_db.entry[i].F;    /* separate point from F ?*/
		double t1 = torc_gettime();
		calculate_statistics(fj, n, data.Num[runinfo.Gen], runinfo.Gen, sel);
		double t2 = torc_gettime();
		printf("init + calc stats : %lf + %lf = %lf seconds\n", t2-t1, t1-t0, t2-t0);
    }

    int newchains = 0;
    for (i = 0; i < n; i++) {
        if (sel[i] != 0) newchains++;
    }

    sort_t *list;
    list = calloc(1, n*sizeof(sort_t));
    for (i = 0; i < n; i++) {
        list[i].idx = i;
        list[i].nsel = sel[i];
        list[i].F = curgen_db.entry[i].F;
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

#if 1   /* UPPER THRESHOLD */
      /* peh:check this */
      /* breaking long chains */
    int initial_newchains = newchains;
    int h_threshold = data.MaxChainLength;    /* peh: configuration file + more balanced lengths */
    for (i = 0; i < initial_newchains; i++) {
        if (list[i].nsel > h_threshold) {
            while (list[i].nsel > h_threshold) {
                list[newchains] = list[i];
                list[newchains].nsel = h_threshold;
                list[i].nsel = list[i].nsel - h_threshold;
                newchains++;
            }
        }
    }

    qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
    printf("Points broken\n");
    for (i = 0; i < n; i++) {
        printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
    }
#endif

#endif

#if 1   /* LOWER THRESHOLD */
    /* single to double step chains */
    int l_threshold = data.MinChainLength;    /* peh: configuration file + more balanced lengths */
    for (i = 0; i < newchains; i++) {
        if ((list[i].nsel > 0)&&(list[i].nsel < l_threshold)) {
            list[i].nsel = l_threshold;
        }
    }

    qsort(list, n, sizeof(sort_t), compar_desc);

#if VERBOSE
    printf("Points advanced\n");
    for (i = 0; i < n; i++) {
        printf("%d: %d %d %f\n", i, list[i].idx, list[i].nsel, list[i].F);
    }
#endif

#endif





    int ldi;    /* leader index*/
    ldi = 0;
    for (i = 0; i < n; i++) {    /* newleader */
        if (list[i].nsel != 0) {
            int idx = list[i].idx;
            for (p = 0; p < data.Nth; p++) {
                leaders[ldi].point[p] = curgen_db.entry[idx].point[p];
            }
            leaders[ldi].F = curgen_db.entry[idx].F;
            leaders[ldi].nsel = list[i].nsel;
            ldi++;
        }
    }

    free(list);

    for (i = 0; i < newchains; i++) leaders[i].queue = -1;    /* rr*/

#if VERBOSE
    printf("Leaders before\n");
    for (i = 0; i < newchains; i++) {
        printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
    }
#endif


#if defined(_USE_TORC_)
    /* cool and greedy partitioning ala Panos-- ;-) */

    int nworkers = torc_num_workers();
    int *workload = (int *)calloc(1, nworkers*sizeof(int));    /* workload[1..workers] = 0*/

    for (i = 0; i < newchains; i++) {
        int least_loader_worker = compute_min_idx_i(workload, nworkers);
        leaders[i].queue = least_loader_worker;
        workload[least_loader_worker] += leaders[i].nsel;
    }

    print_matrix_i((char *)"initial workload", workload, nworkers);
    free(workload);
#endif

#if VERBOSE
    printf("Leaders after\n");
    for (i = 0; i < newchains; i++) {
        printf("%d %d %f %d\n", i, leaders[i].nsel, leaders[i].F, leaders[i].queue);
    }
#endif

    if (1)
    {/*start block*/
        /*    double x[data.Nth][n];*/
        double **x = g_x;
        for (i = 0; i < newchains; i++) {
            for (p = 0; p < data.Nth; p++) {
                x[p][i] = leaders[i].point[p];
            }
        }

        double meanx[data.Nth], stdx[data.Nth];
        for (p = 0; p < data.Nth; p++) {
            meanx[p] = compute_mean(x[p], newchains);
            stdx[p] = compute_std(x[p], newchains, meanx[p]);
        }

        printf("CURGEN DB (LEADER) %d: [nlead=%d]\n", runinfo.Gen, newchains);
		if(data.options.Display){
			print_matrix((char *)"means", meanx, data.Nth);
	        print_matrix((char *)"std", stdx, data.Nth);
		}
    }/*end block*/

    if (data.use_local_cov)
        precompute_chain_covariances(leaders, data.init_mean, data.local_cov, newchains);

    curgen_db.entries = 0;    /* reset curgen db*/
    printf("calculate_statistics: newchains=%d\n", newchains);

    for (i = 0; i < data.Nth; i++) free(g_x[i]);
    free(g_x);

    free(fj);
    free(sel);

    return newchains;
}



