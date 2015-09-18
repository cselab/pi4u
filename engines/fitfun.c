/*
 *  fitfun.c
 *  Pi4U
 *
 *  Created by Lina Kulakova on 16/9/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#define _XOPEN_SOURCE 600
#define _BSD_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <ftw.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>
#include <torc.h>
#include "spawner.c"
#include "gsl_headers.h"

int dbg_display = 1;

char bindir[256];
char datdir[256];
int input_method;		// 0: argv, 1: text file, 2: binary file
int output_method;		// 1: text file, 2: binary file 

void sighandler(int signo)
{
	printf("spawner: SIGNAL %d was received!\n", signo); fflush(0);
}

void fitfuntask(double in_tparam[], int *pdim, double *out_tparam, int winfo[4])
{
	int me = getpid();	/* spanwer_id : worker_id */

	double aij = in_tparam[0];
	double gammadpd = in_tparam[1];
	double v = in_tparam[2];

	char line[1024];
	char args[1024];
	char *largv[64];
	char out_file[1024];

	int rf = fork();
	if (rf < 0) {
		printf("spanwer(%d): fork failed!!!!\n", me); fflush(0);
	}

	if (rf == 0) {
		strcpy(line, "");
		strcat(line, "./test_sphere 1 1 1 -tend=0.002");
		strcpy(args, "");
		sprintf(args, "-aij=%lf -gammadpd=%lf -v=%lf", aij, gammadpd, v);
		strcat(line, args);
	
		parse(line, largv);	/* prepare argv */
	
		sprintf(out_file, "output_v_%.16lf.txt", v);
		int fd = open(out_file, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient
	
		execvp(*largv, largv);
	
		/* TODO: get the Cd */
	}
	waitpid(rf, NULL, 0);
	sync();
}

double fitfun(double *input, int n, void *output, int *info)
{
	double res = 1e8;
	int i, j;
	int me = getpid();	/* spanwer_id : worker_id */
	int rf;
	char line[1024];
	char args[1024];
	char *largv[64];
	char taskname[256];
	double t0, t1;
	double v0 = 1.0;
	double L = 1.0;
	double nu, Re0;
	int winfo[4];
	double in_tparam[3];
	int pdim = 3;

	double aij = input[0];
	double gammadpd = input[1];
	double sigma = input[2];
	double len = input[3];

	double sigma2 = sigma*sigma;
	double len2 = len*len;

	sprintf(taskname, "tmpdir.%d.%d.%d.%d", info[0], info[1], info[2], info[3]);

	if (dbg_display) {
		printf("spanwer(%d): running task %s with params (", me, taskname);

		for (i = 0; i < n-1; i++)
			printf("%.16lf,", input[i]);

		printf("%.16lf)\n", input[i]);

		fflush(0);
	}

	/* create a (temporary) directory */
	mkdir(taskname, S_IRWXU);

	t0 = my_gettime();
	rf = fork();
	if (rf < 0) {
		printf("spanwer(%d): fork failed!!!!\n", me); fflush(0);
	}

	if (rf == 0) {
		chdir(taskname);	/* enter to the new directory */

		/* 1. PREPROCESSING PHASE - APPLICATION SPECIFIC */
		/* copy some application specific input files*/
		copy_file("/users/lina/UQ/pi4u.git/engines", "drag_vs_re_data.txt");
		copy_file("/users/lina/UQ/pi4u.git/engines/fitfun_code/mpi-dpd", "test");
#if 0
		copy_file("/users/lina/UQ/pi4u.git/engines/fitfun_code/mpi-dpd", "test_sphere");
#endif

		/* 1a. Write parameters to the file */
		FILE *finp = fopen("params.dat", "w");
		for (i = 0; i < n; i++)
			fprintf(finp, "%.16lf\n", input[i]);
		fclose(finp);

		/* 1b. Compute viscosity */
		strcpy(line, "");
		strcat(line, "./test 1 1 1 -tend=0.002"); /* TODO: check command line options */
		strcpy(args, "");
		sprintf(args, " -aij=%lf -gammadpd=%lf -v=%lf", aij, gammadpd, v0);
		strcat(line, args);

#if 1
		printf("Running command: %s\n", line);
		parse(line, largv);	/* prepare argv */
		printf("execvp command:");
		for(i=0; ; ++i)
		{
			if(largv[i] != '\0')
				printf(" %s", largv[i]);
			else
			{
				printf("\n", largv[i]);
				break;
			}
		}
#endif

		int fd = open("output_visc.txt", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient

/*		char *argm[] = {"ls", "-la", 0}; 
		execvp(argm[0], argm); */
		execvp(*largv, largv);

		/* TODO: get the viscosity */
		nu = 1.0;
		Re0 = v0*L/nu;
	}
	waitpid(rf, NULL, 0);
	sync();

#if 0
	/* 2. EXECUTION OF THE FUNCTION EVALUATION (SIMULATION) SOFTWARE */
	/* 2a. Read data from file */

	FILE * fdata = fopen("drag_vs_re_data.txt", "r");
	int ndata;
	fscanf(fdata, "%d\n", &ndata);
	double * Re = malloc(ndata*sizeof(double));
	double * Cd_exp = malloc(ndata*sizeof(double));
	for(i=0; i<ndata; ++i)
		fscanf(fdata, "%lf %lf\n", &Re[i], &Cd_exp[i]);
	fclose(fdata);

	/* 2b. Compute drag vs Reynolds */
	double * Cd_sim = malloc(ndata*sizeof(double));
	for(i=0; i<ndata; ++i)
	{
		/* Run the simulation with current Re */
		winfo[0] = info[0]; /* generation id */
		winfo[1] = info[1]; /* chain id */
		winfo[2] = info[2]; /* task in chain id */
		winfo[3] = i;       /* subtask id */

		in_tparam[0] = aij;
		in_tparam[1] = gammadpd;
		in_tparam[2] = Re[i]/(Re0*v0); /* v */
		torc_create(-1, fitfuntask, 4,
			pdim,	MPI_DOUBLE, CALL_BY_COP,
			1,		MPI_INT,	CALL_BY_COP,
			1,		MPI_DOUBLE, CALL_BY_RES,
			4,		MPI_INT,	CALL_BY_COP,
			in_tparam, &pdim, &Cd_sim[i], winfo);
	}
	/* wait for all chain tasks to finish */
#ifdef _STEALING_
		torc_enable_stealing();
#endif
		torc_waitall();
#ifdef _STEALING_
		torc_disable_stealing();
#endif

	t1 = my_gettime();
		
	if (dbg_display)
		printf("spanwer(%d): simcode_time=%lf secs\n", me, t1-t0);fflush(0);

	/* 3. POSTPROCESSING PHASE - APPLICATION SPECIFIC */
	/* write results to file */
	FILE * fres = fopen("drag_vs_re_result.txt", "r");
	fprintf(fres, "%d\n", ndata);
	for(i=0; i<ndata; ++i)
		fprintf(fres, "%lf %lf\n", Re[i], Cd_sim[i]);
	fclose(fres);

	/* compute fitness */
	gsl_vector * diff = gsl_vector_calloc(ndata);
	gsl_matrix * cov = gsl_matrix_calloc(ndata, ndata);
	gsl_matrix * cov_inv = gsl_matrix_calloc(ndata, ndata);
	gsl_permutation * perm = gsl_permutation_alloc(ndata);
	gsl_vector * temp = gsl_vector_calloc(ndata);
	int signum;
	double res_cov_part;

	for(i=0; i<ndata; ++i)
	{
		gsl_vector_set(diff, i, Cd_exp[i]-Cd_sim[i]);
		gsl_matrix_set(cov, i, i, sigma2);
		for(j=i+1; j<ndata; ++j)
		{
			double curr = sigma2*exp(-0.5*pow(Re[i]-Re[j],2)/len2);
			gsl_matrix_set(cov, i, j, curr);
			gsl_matrix_set(cov, j, i, curr);
		}
	}

	gsl_linalg_LU_decomp(cov, perm, &signum); /* now cov contains its LU decomposition */
	gsl_linalg_LU_invert(cov, perm, cov_inv); /* cov_inv contains inverse of cov */

	double log_det_cov = gsl_linalg_LU_lndet(cov);
	gsl_blas_dgemv(CblasNoTrans, 1.0, cov_inv, diff, 0.0, temp);
	gsl_blas_ddot(diff, temp, &res_cov_part);

	res = -0.5*( ndata*log(2*M_PI) + log_det_cov + res_cov_part);

	gsl_vector_free(diff);
	gsl_matrix_free(cov);
	gsl_matrix_free(cov_inv);
	gsl_permutation_free(perm);
	gsl_vector_free(temp);

	/* remove the temporary directory */
	if (0)
		rmrf(taskname);

	free(Re);
	free(Cd_exp);
	free(Cd_sim);

#endif

	return res;
}

