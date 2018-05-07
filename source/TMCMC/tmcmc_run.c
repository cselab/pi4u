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

#include "tmcmc_db.h"
#include "tmcmc_aux.h"
#include "tmcmc_stats.h"
#include "tmcmc_engine.h"
#include "fitfun.h"


#include "../priors/priors.h"
#include "../priors/myrand.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>


#if defined(_USE_TORC_)
	#include <mpi.h>
	#ifdef __cplusplus
		extern "C"
		{
	#endif
	#include <torc.h>
	#ifdef __cplusplus
		}
	#endif
#else

	#include <pthread.h>

	#if defined(_USE_OPENMP_)
		#include <omp.h>
	#endif

#endif



extern Density *priors;




int main(int argc, char *argv[])
{

    int i;
    double t0, gt0, gt1;
    int winfo[4];
    int nchains; 

	#if defined(_USE_TORC_)
    	torc_register_task((void *)initchaintask);
		torc_register_task((void *)chaintask);
		torc_register_task((void *)torc_update_full_db_task);
		torc_register_task((void *)torc_update_curgen_db_task);
		torc_register_task((void *)torc_update_curres_db_task);
		torc_register_task((void *)reset_nfc_task);
		torc_register_task((void *)get_nfc_task);
		torc_register_task((void *)taskfun);
		torc_register_task((void *)call_gsl_rand_init);
		torc_register_task((void *)call_print_matrix_2d);
		torc_register_task((void *)call_update_gdata);
	#endif

    data_init();

	int Npr;
	read_priors( "priors.par", &priors, &Npr );
	print_priors( priors, Npr);
	//XXX: check that Nps==data.Nth




    char str[12];
    sprintf(str, "%d", data.Nth);
    fitfun_initialize(str);

	#if defined(_USE_TORC_)
		torc_init(argc, argv, MODE_MS);
	#endif
	
	#if defined(_AFFINITY_)
		spmd_setaffinity();
	#endif



    spmd_gsl_rand_init();


    curgen_db.entries = 0; 

	//-----  Restart: load data -------------------------
    int goto_next = 0;
    int res;
    if (data.restart) {
    	res = load_runinfo();
    	if (res == 0) {
        	load_curgen_db(runinfo.Gen);
        	nchains = data.Num[0];
        	printf("nchains = %d\n", nchains);
        	gt0 = t0 = torc_gettime();
        	goto_next = 1;
    	}
    }
	//---------------------------------------------------




    gt0 = t0 = torc_gettime();

    nchains = data.Num[0];
    double out_tparam[data.PopSize];    /* nchains*/

    
	// sample from the prior
	if (goto_next == 0){
       

		#if defined(_USE_OPENMP_)
			#pragma omp parallel
			{
				printf("Hello from thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
				#pragma omp for
			//	{
		#endif


        for (i = 0; i < nchains; i++){
            winfo[0] = runinfo.Gen;
            winfo[1] = i;
            winfo[2] = -1;
            winfo[3] = -1;

            double in_tparam[data.Nth]; // here store the samples
			for (int d = 0; d < data.Nth; d++)
				in_tparam[d] = eval_random( priors[d] );

			#if defined(_USE_TORC_)
                torc_create(-1, (void (*)())initchaintask, 4,
                        data.Nth, MPI_DOUBLE, CALL_BY_COP,
                        1, MPI_INT, CALL_BY_COP,
                        1, MPI_DOUBLE, CALL_BY_RES,
                        4, MPI_INT, CALL_BY_COP,
                        in_tparam, &data.Nth, &out_tparam[i], winfo);
			#else
				
				
				#if defined(_USE_OPENMP_)
					//#pragma omp task shared(data, out_tparam) firstprivate(i, in_tparam, winfo)
              		#pragma omp task firstprivate(i, winfo, in_tparam) shared(data, out_tparam)
				#endif
				{
            		initchaintask(in_tparam, &data.Nth, &out_tparam[i], winfo);
                }
			
			#endif

        }

		#if defined(_USE_OPENMP_)
			//} // single
			} 	// parallel
		#endif


		#if defined(_USE_TORC_)
        	if (data.stealing)
	    		torc_enable_stealing();

        	torc_waitall();

        	if (data.stealing)
            	torc_disable_stealing();
		#endif


        gt1 = torc_gettime();
        int g_nfeval = get_nfc();
        printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", 
																							runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
        reset_nfc();
	
        if( data.icdump ) dump_curgen_db(runinfo.Gen);
        if( data.ifdump ) dump_full_db(runinfo.Gen);

		save_runinfo();
		if (data.restart)	check_for_exit();

        if (data.MaxStages == 1) goto end;

    }





    static cgdbp_t *leaders; /*[MAXCHAINS];*/
    leaders = (cgdbp_t *) calloc( 1, data.PopSize*sizeof(cgdbp_t) );
    for (i = 0; i < data.PopSize; i++){
        leaders[i].point = (double *) calloc(1, data.Nth*sizeof(double));
    }

    curres_db.entries = 0;
    nchains = prepare_newgen(nchains, leaders);    /* calculate statistics */

    spmd_update_gdata();
	if(data.options.Display)
		call_print_matrix_2d();

    /* this can be moved above */
    if (runinfo.p[runinfo.Gen] == 1) {
        printf("p == 1 from previous run, nothing more to do\n");
        goto end;
    }

	printf("----------------------------------------------------------------\n");

    runinfo.Gen++;
    for (; runinfo.Gen < data.MaxStages; runinfo.Gen++){

        int nsteps;
        gt0 = torc_gettime();


		#if defined(_USE_OPENMP_)
			#pragma omp parallel
			{
			#pragma omp single nowait
			{
		#endif
        
		int winfo[4];
        double in_tparam[data.Nth];
        double init_mean[data.Nth];
        double chain_cov[data.Nth*data.Nth];

        for (i = 0; i < nchains; i++) {
            winfo[0] = runinfo.Gen;
            winfo[1] = i;
            winfo[2] = -1;    /* not used */
            winfo[3] = -1;    /* not used */

            int p;
            for (p = 0; p < data.Nth; p++)
                in_tparam[p] = leaders[i].point[p];
            
			nsteps = leaders[i].nsel;

            if (data.use_local_cov){
                for (p = 0; p < data.Nth*data.Nth; ++p)
                    chain_cov[p] = data.local_cov[i][p];

                for (p = 0; p < data.Nth; ++p) {
                    if (data.use_proposal_cma)
                        init_mean[p] = data.init_mean[i][p];
                    else
                        init_mean[p] = leaders[i].point[p];
                }
            }
            else{
                for (p = 0; p < data.Nth; ++p)
                    for (int j = 0; j < data.Nth; ++j)
                        chain_cov[p*data.Nth+j]= data.bbeta*runinfo.SS[p][j];

                for (p = 0; p < data.Nth; ++p)
                    init_mean[p] = in_tparam[p];
            }

            out_tparam[i] = leaders[i].F;    /* loglik_leader...*/

			#if defined(_USE_TORC_)
            	torc_create(leaders[i].queue, (void (*)())chaintask, 7,
                    				data.Nth, MPI_DOUBLE, CALL_BY_COP,
                    				1, MPI_INT, CALL_BY_COP,
                    				1, MPI_INT, CALL_BY_COP,
                    				1, MPI_DOUBLE, CALL_BY_REF,
                    				4, MPI_INT, CALL_BY_COP,
                    				data.Nth, MPI_DOUBLE, CALL_BY_COP,
                    				data.Nth*data.Nth, MPI_DOUBLE, CALL_BY_COP,
                    				in_tparam, &data.Nth, &nsteps, &out_tparam[i], winfo,
                    				init_mean, chain_cov);
			#else

				#if defined(_USE_OPENMP_)
	    			#pragma omp task shared(data, out_tparam) firstprivate(i, nsteps, in_tparam, winfo, init_mean, chain_cov)
				#endif
            	chaintask( in_tparam, &data.Nth, &nsteps, &out_tparam[i], winfo, init_mean, chain_cov );

			#endif
        }

        // wait for all chain tasks to finish
		#if defined(_USE_OPENMP_)
			} // single
			} // parallel
		#endif


		#if defined(_USE_TORC_)
			if (data.stealing)
            	torc_enable_stealing();

        	torc_waitall();
        	if (data.stealing)
            	torc_disable_stealing();
		#endif

		

        gt1 = torc_gettime();
        int g_nfeval = get_nfc();
        printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", 
																	runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
        reset_nfc();
		


        if (data.icdump) dump_curgen_db(runinfo.Gen);
        if (data.ifdump) dump_full_db(runinfo.Gen);

        // save here
		save_runinfo();
        if (data.restart) {
            check_for_exit();
        }

        curres_db.entries = 0;
        nchains = prepare_newgen(nchains, leaders);    // calculate statistics

        spmd_update_gdata();

		if(data.options.Display)
			call_print_matrix_2d();

		printf("Acceptance rate   :  %lf \n",runinfo.acceptance[runinfo.Gen]) ;
		printf("Annealing exponent:  %lf \n",runinfo.p[runinfo.Gen]) ;
		printf("----------------------------------------------------------------\n");


        if (runinfo.p[runinfo.Gen] == 1) break;

        if (runinfo.Gen+1 == data.MaxStages) break;

    } // loop over generations

    if (data.MaxStages == 1) runinfo.Gen = 0;	// small correction for this extreme case


	//FIXME move this prints into a file
    print_matrix((char *)"runinfo.p", runinfo.p, runinfo.Gen+1);
    print_matrix((char *)"runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
    print_matrix_i((char *)"runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
    print_matrix((char *)"runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
    print_matrix((char *)"runinfo.logselection", runinfo.logselection, runinfo.Gen+1);

    double logEvidence[1];
    logEvidence[0] = compute_sum(runinfo.logselection, runinfo.Gen+1);
    print_matrix((char *)"logEvidence", logEvidence, 1);


	FILE *fp;
	fp = fopen("log_evidence.txt", "w");
	fprintf(fp, "%lf\n", logEvidence[0]);
	fclose(fp);


    save_runinfo();




end:
    /* making a copy of last curgen_db file */
    if (data.icdump){
        printf("lastgen = %d\n", runinfo.Gen);
        char cmd[256];
        sprintf(cmd, "cp curgen_db_%03d.txt final.txt", runinfo.Gen);
        system(cmd);
    }


    /* shutdown */
    fitfun_finalize();

    printf("total function calls = %d\n", get_tfc());
	#if defined(_USE_TORC_)
    	torc_finalize();
	#endif
    
	
		
	return 0;
}
