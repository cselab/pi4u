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


#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


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

    gt0 = t0 = torc_gettime();

    nchains = data.Num[0];
    double out_tparam[data.PopSize];    /* nchains*/

    if (goto_next == 0)
    {
        FILE *init_fp = NULL;
        if (data.load_from_file == 1) {
            init_fp = fopen("init_db.txt", "r");    /* peh: parametrize this */
            if (init_fp == NULL) {
                printf("init_db.txt file not found!\n");
                exit(1);
            }
        }

		#if defined(_USE_OPENMP_)
		#pragma omp parallel
			{
			printf("Hello from thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
		#pragma omp for
		//	{
		#endif

        for (i = 0; i < nchains; i++) {
            winfo[0] = runinfo.Gen;
            winfo[1] = i;
            winfo[2] = -1;
            winfo[3] = -1;

            double in_tparam[data.Nth];

            if (data.prior_type == 0)    /* uniform */
            {
                /* uniform */
                int d;
                for (d = 0; d < data.Nth; d++) {
                    in_tparam[d] = uniformrand(0,1);
                    in_tparam[d] *= (data.upperbound[d]-data.lowerbound[d]);
                    in_tparam[d] += data.lowerbound[d];
                }
            }
            else if (data.prior_type == 1)    /* gaussian */
            {
                mvnrnd(data.prior_mu, data.prior_sigma, in_tparam, data.Nth);
            }
            else if (data.prior_type == 3)    /* composite */
            {
                int d;
                for (d = 0; d < data.Nth; d++) {
                    if (data.compositeprior_distr[d] == 0) {
                        in_tparam[d] = uniformrand(0,1);
                        in_tparam[d] *= (data.upperbound[d]-data.lowerbound[d]);
                        in_tparam[d] += data.lowerbound[d];
                    }
                    else if (data.compositeprior_distr[d] == 1) {
                        mvnrnd(&data.prior_mu[d], &data.prior_sigma[d], &in_tparam[d], 1);
                    }
                    else if (data.compositeprior_distr[d] == 2) {
                      	in_tparam[d] = truncated_normal_rand (data.prior_mu[d], data.prior_sigma[d], data.lowerbound[d], data.upperbound[d]);
                    }
                }
            }

            if (data.load_from_file == 1)    /* override the computed point and read it from the file */
            {
              	printf("reading from init_db.txt\n");
                int j;
                for (j = 0; j < data.Nth; j++) fscanf(init_fp, "%lf", &in_tparam[j]);
                fscanf(init_fp, "%lf", &out_tparam[i]);
                double prior_val;
                fscanf(init_fp, "%lf", &prior_val);

                /*torc_update_curgen_db(in_tparam, out_tparam[i]);*/    /* peh - eval or not */
                /*if (data.ifdump) torc_update_full_db(in_tparam, out_tparam[i], NULL, 0, 0);*/
            }

            if (data.prior_type <= 3)    /* peh: file without or with function evaluations? */
            {
#if defined(_USE_TORC_)
                torc_create(-1, (void (*)())initchaintask, 4,
                        data.Nth, MPI_DOUBLE, CALL_BY_COP,
                        1, MPI_INT, CALL_BY_COP,
                        1, MPI_DOUBLE, CALL_BY_RES,
                        4, MPI_INT, CALL_BY_COP,
                        in_tparam, &data.Nth, &out_tparam[i], winfo);
#else

//#if defined(_USE_OPENMP_)
//#endif

//                #pragma omp task shared(data, out_tparam) firstprivate(i, in_tparam, winfo)
//                {
//                initchaintask(
//                        in_tparam, &data.Nth, &out_tparam[i], winfo);
//                }

                #pragma omp task firstprivate(i, winfo, in_tparam) shared(data, out_tparam)
                {
                  //usleep(10*1000);
                  //printf("hello from task %d on worker %d\n", i, omp_get_thread_num());
                  initchaintask(in_tparam, &data.Nth, &out_tparam[i], winfo);
                }

#endif
            }
        }

#if defined(_USE_OPENMP_)
	//} // single
	} // parallel
#endif


#if defined(_USE_TORC_)
        if (data.stealing)
	    torc_enable_stealing();

        torc_waitall();

        if (data.stealing)
            torc_disable_stealing();
#endif

        if (data.load_from_file == 1) {
            fclose(init_fp);
        }

        gt1 = torc_gettime();
        int g_nfeval = get_nfc();
        printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
        reset_nfc();

        if (data.icdump) dump_curgen_db(runinfo.Gen);
        if (data.ifdump) dump_full_db(runinfo.Gen);

        /* save here */
		save_runinfo();
		if (data.restart) {
            check_for_exit();
        }

        if (data.MaxStages == 1) goto end;

    }
    ;

    static cgdbp_t *leaders; /*[MAXCHAINS];*/
    leaders = (cgdbp_t *)calloc(1, data.PopSize*sizeof(cgdbp_t));
    for (i = 0; i < data.PopSize; i++) {
        leaders[i].point = (double *)calloc(1, data.Nth*sizeof(double));
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

        /* process current generation, compute probs, find new chains */
        /*leader[i]: { point[data.Nth], F, nsteps}*/

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

            if (data.use_local_cov) {
                for (p = 0; p < data.Nth*data.Nth; ++p)
                    chain_cov[p] = data.local_cov[i][p];

                for (p = 0; p < data.Nth; ++p) {
                    if (data.use_proposal_cma)
                        init_mean[p] = data.init_mean[i][p];
                    else
                        init_mean[p] = leaders[i].point[p];
                }
            }
            else {
                int j;
                for (p = 0; p < data.Nth; ++p)
                    for (j = 0; j < data.Nth; ++j)
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
            chaintask(
                    in_tparam, &data.Nth, &nsteps, &out_tparam[i], winfo,
                    init_mean, chain_cov);
#endif
        }

        /* wait for all chain tasks to finish */

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
        printf("server: Generation %d: total elapsed time = %lf secs, generation elapsed time = %lf secs for function calls = %d\n", runinfo.Gen, gt1-t0, gt1-gt0, g_nfeval);
        reset_nfc();

        if (data.icdump) dump_curgen_db(runinfo.Gen);
        if (data.ifdump) dump_full_db(runinfo.Gen);

        /* save here*/
		save_runinfo();
        if (data.restart) {
            check_for_exit();
        }

        curres_db.entries = 0;
        nchains = prepare_newgen(nchains, leaders);    /* calculate statistics*/

        spmd_update_gdata();
        /*spmd_print_matrix_2d();*/
		if(data.options.Display)
			call_print_matrix_2d();

		printf("Acceptance rate   :  %lf \n",runinfo.acceptance[runinfo.Gen]) ;
		printf("Annealing exponent:  %lf \n",runinfo.p[runinfo.Gen]) ;
		printf("----------------------------------------------------------------\n");


#if 0
        printf("=================\n");
        print_matrix((char *)"runinfo.p", runinfo.p, runinfo.Gen+1);
        print_matrix((char *)"runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
        print_matrix_i((char *)"runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
        print_matrix((char *)"runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
        print_matrix((char *)"runinfo.logselection", runinfo.logselection, runinfo.Gen+1);
        printf("=================\n");
#endif

        if (runinfo.p[runinfo.Gen] == 1) {
            break;
        }
        if (runinfo.Gen+1 == data.MaxStages) {
            break;
        }
    }

    if (data.MaxStages == 1) runinfo.Gen = 0;	// small correction for this extreme case

    print_matrix((char *)"runinfo.p", runinfo.p, runinfo.Gen+1);
    print_matrix((char *)"runinfo.CoefVar", runinfo.CoefVar, runinfo.Gen+1);
    print_matrix_i((char *)"runinfo.currentuniques", runinfo.currentuniques, runinfo.Gen+1);
    print_matrix((char *)"runinfo.acceptance", runinfo.acceptance, runinfo.Gen+1);
    print_matrix((char *)"runinfo.logselection", runinfo.logselection, runinfo.Gen+1);

    double logEvidence[1];
    logEvidence[0] = compute_sum(runinfo.logselection, runinfo.Gen+1);
    print_matrix((char *)"logEvidence", logEvidence, 1);

    /* peh:check -- inner tmcmc */
    {
        FILE *fp;
        fp = fopen("fitness.txt", "w");
        fprintf(fp, "%lf\n", logEvidence[0]);
        fclose(fp);
    }

    // print_matrix_2d((char *)"runinfo.SS", runinfo.SS, data.Nth, data.Nth);

    // for (i = 0; i < runinfo.Gen+1; i++) {
    //     char title[64];
    //     sprintf(title, "runinfo.meantheta(%d)", i);
    //     /*print_matrix((char *)"runinfo.meantheta", runinfo.meantheta[i], data.Nth);*/
    //     print_matrix((char *)title, runinfo.meantheta[i], data.Nth);
    // }

    /* last save here - do we need this? what happens if I restart the program with this saved data*/
    // if (data.restart)
    save_runinfo();

end:
    /* making a copy of last curgen_db file */
    if (data.icdump)
    {
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
