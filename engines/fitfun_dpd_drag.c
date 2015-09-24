/*
 *  fitfun.c
 *  Pi4U
 *
 *  Created by Lina Kulakova on 16/9/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#define _XOPEN_SOURCE 700
#define _BSD_SOURCE 1
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

void fitfuntask(double in_tparam[], int *pdim, double *out_tparam, char taskname[], int winfo[4])
{
	int me = getpid();	/* spanwer_id : worker_id */

	double aij = in_tparam[0];
	double gammadpd = in_tparam[1];
	double v = in_tparam[2];

	char line[1024];
	char args[1024];
	char *largv[64];
	char out_file[1024];

	printf("Generation: %d, chain: %d, task: %d\n", winfo[0], winfo[1], winfo[2]);

	int rf = fork();
	if (rf < 0) {
		printf("spanwer(%d): fork failed!!!!\n", me); fflush(0);
	}

	if (rf == 0) {
		chdir(taskname);	/* enter to the new directory */

		strcpy(line, "");
		strcat(line, "./test 1 1 1 -tend=0.002"); /* TODO: fix executable */
		strcpy(args, "");
		sprintf(args, " -aij=%lf -gammadpd=%lf -v=%lf", aij, gammadpd, v);
		strcat(line, args);

		sprintf(out_file, "output_v_%.16lf.txt", v);
		int fd = open(out_file, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient
#if 1
		parse(line, largv);	/* prepare argv */
		execvp(*largv, largv);
#else
		printf("Running command for Re vs drag %s\n", line);
#endif
	}
	waitpid(rf, NULL, 0);
	sync();
	
	/* TODO: get the Cd */
	*out_tparam = 3e6;
}

double fitfun(double *input, int n, void *output, int *info)
{
	const char * basedir = "/users/lina/UQ/pi4u.git/engines";
	const char *data_filename = "drag_vs_re_data.txt";
	const char *res_filename = "drag_vs_re_result.txt";
	double res = 1e8;
	int i, j;
	int me = getpid();	/* spanwer_id : worker_id */
	int rf;
	char buf[1024];
	char line[1024];
	char args[1024];
	char *largv[64];
	char taskname[256];
	double t0, t1;
	double v0 = 1.0;
	double L = 1.0; /* TODO: check cylinder diameter */
	double nu = 1.0, Re0 = 1.0;
	int winfo[4];
	const int pdim = 3;
	double in_tparam[pdim];

	double aij = input[0];
	double gammadpd = input[1];
	double sigma = input[2];
	double len = input[3];

	double sigma2 = sigma*sigma;
	double inv_len2 = 1.0/(len*len);

	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

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
	/* 1. PREPROCESSING PHASE - APPLICATION SPECIFIC */
	/* copy some application specific input files*/
	copy_file(basedir, data_filename);

	rf = fork();
	if (rf < 0) {
		printf("spanwer(%d): fork failed!!!!\n", me); fflush(0);
	}

	if (rf == 0) {
		chdir(taskname);	/* enter to the new directory */

		strcpy(buf, "");
		sprintf(buf, "%s/fitfun_code/mpi-dpd", basedir);
		copy_file(buf, "test");
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
		int fd = open("output_visc.txt", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient
#if 1
		parse(line, largv);	/* prepare argv */
		execvp(*largv, largv);
#else
		printf("Running command for viscosity: %s\n", line);
#endif
	}
	waitpid(rf, NULL, 0);
	sync();

	/* TODO: get the viscosity */
	nu = 1.0;
	Re0 = v0*L/nu;

	/* 2. EXECUTION OF THE FUNCTION EVALUATION (SIMULATION) SOFTWARE */
	/* 2a. Read data from file */

	strcpy(buf, "");
	sprintf(buf, "%s/%s", workdir, data_filename);
	FILE * fdata = fopen(buf, "r");
	if(!fdata)
	{
		printf("Cannot open file %s\n", buf); fflush(0);
		return res;
	}
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
		winfo[2] = i;		/* task id */
		winfo[3] = 0;       /* unused */

		in_tparam[0] = aij;
		in_tparam[1] = gammadpd;
		in_tparam[2] = Re[i]/(Re0*v0); /* v */
		torc_create(-1, fitfuntask, 5,
			pdim,				MPI_DOUBLE, CALL_BY_COP,
			1,					MPI_INT,	CALL_BY_COP,
			1,					MPI_DOUBLE, CALL_BY_RES,
			strlen(taskname),	MPI_INT,	CALL_BY_COP, /* TODO: check with Panos if MPI_CHAR is supported*/
			4,					MPI_INT,	CALL_BY_COP,
			in_tparam, &pdim, &Cd_sim[i], taskname, winfo);
	}
#ifdef _STEALING_
	torc_enable_stealing();
#endif
	/* wait for all chain tasks to finish */
	torc_waitall();
#ifdef _STEALING_
	torc_disable_stealing();
#endif

	t1 = my_gettime();
		
	if (dbg_display)
		printf("spanwer(%d): simcode_time=%lf secs\n", me, t1-t0);fflush(0);

	/* 3. POSTPROCESSING PHASE - APPLICATION SPECIFIC */
	/* write results to file */
	strcpy(buf, "");
	sprintf(buf, "%s/%s", taskname, res_filename);
	FILE * fres = fopen(buf, "w");
	if(!fres)
		printf("Cannot open file %s\n", buf); fflush(0);
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
	double log_cov_part;

	for(i=0; i<ndata; ++i)
	{
		gsl_vector_set(diff, i, Cd_exp[i]-Cd_sim[i]);
		gsl_matrix_set(cov, i, i, sigma2);
		for(j=i+1; j<ndata; ++j)
		{
			double curr = sigma2*exp(-0.5*pow(Re[i]-Re[j],2)*inv_len2);
			gsl_matrix_set(cov, i, j, curr);
			gsl_matrix_set(cov, j, i, curr);
		}
	}

	gsl_linalg_LU_decomp(cov, perm, &signum); /* now cov contains its LU decomposition */
	gsl_linalg_LU_invert(cov, perm, cov_inv); /* cov_inv contains inverse of cov */

	double log_det_cov = gsl_linalg_LU_lndet(cov);
	gsl_blas_dgemv(CblasNoTrans, 1.0, cov_inv, diff, 0.0, temp);
	gsl_blas_ddot(diff, temp, &log_cov_part);

	res = -0.5*( ndata*log(2*M_PI) + log_det_cov + log_cov_part);

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

	return res;
}

