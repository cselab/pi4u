/*
 *  fitfun.c
 *  Pi4U
 *
 *  Created by Lina Kulakova on 16/9/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

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
#include "gsl_headers.h"
#include "spawner.h"

int dbg_display = 1;

char bindir[256];
char datdir[256];
int input_method;		// 0: argv, 1: text file, 2: binary file
int output_method;		// 1: text file, 2: binary file 

void sighandler(int signo)
{
	printf("spawner: SIGNAL %d was received!\n", signo); fflush(0);
}

void fitfuntask(double in_tparam[], int *pdim, double *out_tparam1, double *out_tparam2,
	char taskname[], char workdir[], int winfo[4])
{
	int me = getpid();	/* spanwer_id : worker_id */

	double aij		= in_tparam[0];
	double gammadpd = in_tparam[1];
	double lmax		= in_tparam[2];
	double p		= in_tparam[3];
	double kb		= in_tparam[4];
	double ka		= in_tparam[5];
	double kv		= in_tparam[6];
	double gammaC	= in_tparam[7];
	double force	= in_tparam[8];

	char line[1024];
	char args[1024];
	char *largv[64];
	char buf[1024];

	char out_file[1024];
	sprintf(out_file, "%s/%s/output_force_%.16lf.txt", workdir, taskname, force);

	int rf = fork();
	if (rf < 0) {
		printf("spanwer(%d): fork failed!!!!\n", me); fflush(0);
	}

	if (rf == 0) {
		chdir(taskname);	/* enter to the new directory */

		strcpy(line, "");
		strcat(line, "./test 1 1 1 -tend=0.002"); /* TODO: fix executable */
		strcpy(args, "");
		sprintf(args, " -aij=%lf -gammadpd=%lf -lmax=%lf -p=%lf -kb=%lf -ka=%lf -kv=%lf -gammaC=%lf -force=%lf", aij, gammadpd, lmax, p, kb, ka, kv, gammaC, force);
		strcat(line, args);

		int fd = open(out_file, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient
#if 1
		parse(line, largv);	/* prepare argv */
		execvp(*largv, largv);
#else
		printf("Running command for diam vs force %s\n", line);
#endif
	}
	waitpid(rf, NULL, 0);
	sync();

	chdir(workdir);
	
	/* get the AD, TD */
	double diam[2][1024];
	int count = 0;

	FILE * fp = fopen(out_file, "r");
	if(!fp)
		printf("Cannot open file %s\n", out_file);
	else
	{
		while(fgets(buf, 1024, fp) != NULL)
		{
			if((strstr(buf, "diameters: ")) != NULL) /* TODO: check which string should be looking for */
			{
				sscanf(buf, "%*s %lf %lf\n", &diam[0][count], &diam[1][count]);
				printf("A match found with diameters: %lf %lf\n", diam[0][count], diam[1][count]);
				count++;
			}
		}
		fclose(fp);
	}

	/* average over the last 10% of the data */
	if(count>0)
	{
		int i, num = 0;
		*out_tparam1 = 0;
		*out_tparam2 = 0;
		for(i=floor(0.9*count); i<count; ++i)
		{
			*out_tparam1 += diam[0][i];
			*out_tparam2 += diam[1][i];
			num++;
		}
		*out_tparam1 /= num;
		*out_tparam2 /= num;
	}
	else
	{
		*out_tparam1 = 4e6;
		*out_tparam2 = 2e6;
	}
}

double fitfun(double *input, int n, void *output, int *info)
{
	const char * basedir = "/users/lina/UQ/pi4u.git/engines";
	const char *data_filename = "diam_vs_force_data.txt";
	const char *res_filename = "diam_vs_force_result.txt";
	double res = 1e8;
	int i, j;
	int me = getpid();	/* spanwer_id : worker_id */
	char buf[1024];
	char taskname[256];
	double t0, t1;
	int winfo[4];
	const int pdim = 9;
	double in_tparam[pdim];

#if 0
	double aij		= input[0];
	double gammadpd = input[1];
	double lmax		= input[2];
	double p		= input[3];
	double kb		= input[4];
	double ka		= input[5];
	double kv		= input[6];
	double gammaC	= input[7];
#endif
	double sigma	= input[8];
	double len		= input[9];

	double sigma2 = sigma*sigma;
	double mult_len = 0.5/(len*len);

	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	sprintf(taskname, "tmpdir.%d.%d.%d.%d", info[0], info[1], info[2], info[3]);

	if (dbg_display) {
		printf("spanwer(%d): running task %s with params (", me, taskname);

		for (i = 0; i < n-1; i++)
			printf("%.16lf, ", input[i]);

		printf("%.16lf)\n", input[i]);

		fflush(0);
	}

	/* create a (temporary) directory */
	mkdir(taskname, S_IRWXU);

	t0 = my_gettime();
	/* 1. PREPROCESSING PHASE - APPLICATION SPECIFIC */
	/* 1a. Copy some application specific input files*/
	copy_file(basedir, data_filename);

	chdir(taskname);	/* enter to the new directory */
	strcpy(buf, "");
	sprintf(buf, "%s/fitfun_code/mpi-dpd", basedir);
	copy_file(buf, "test"); /* TODO: change executable name */

	/* 1b. Write parameters to the file */
	FILE *finp = fopen("params.dat", "w");
	for (i = 0; i < n; i++)
		fprintf(finp, "%.16lf\n", input[i]);
	fclose(finp);
	chdir(workdir);

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
	double * force = malloc(ndata*sizeof(double));
	double * AD_exp = malloc(ndata*sizeof(double));
	double * TD_exp = malloc(ndata*sizeof(double));
	for(i=0; i<ndata; ++i)
		fscanf(fdata, "%lf %lf %lf\n", &force[i], &AD_exp[i], &TD_exp[i]);
	fclose(fdata);

	/* 2b. Compute stretching vs force */
	double * AD_sim = malloc(ndata*sizeof(double));
	double * TD_sim = malloc(ndata*sizeof(double));
	for(i=0; i<ndata; ++i)
	{
		/* Run the simulation with current force */
		winfo[0] = info[0]; /* generation id */
		winfo[1] = info[1]; /* chain id */
		winfo[2] = i;		/* task id */
		winfo[3] = 0;       /* unused */

		in_tparam[0] = input[0];
		in_tparam[1] = input[1];
		in_tparam[2] = input[2];
		in_tparam[3] = input[3];
		in_tparam[4] = input[4];
		in_tparam[5] = input[5];
		in_tparam[6] = input[6];
		in_tparam[7] = input[7];
		in_tparam[8] = force[i];
		torc_create(-1, fitfuntask, 7,
			pdim,				MPI_DOUBLE, CALL_BY_COP,
			1,					MPI_INT,	CALL_BY_COP,
			1,					MPI_DOUBLE, CALL_BY_RES,
			1,					MPI_DOUBLE, CALL_BY_RES,
			strlen(taskname),	MPI_INT,	CALL_BY_COP, /* TODO: check with Panos if MPI_CHAR is supported*/
			strlen(workdir),	MPI_INT,	CALL_BY_COP, /* TODO: check with Panos if MPI_CHAR is supported*/
			4,					MPI_INT,	CALL_BY_COP,
			in_tparam, &pdim, &AD_sim[i], &TD_sim[i], taskname, workdir, winfo);
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
	sprintf(buf, "%s/%s/%s", workdir, taskname, res_filename);
	FILE * fres = fopen(buf, "w");
	if(!fres)
		printf("Cannot open file %s\n", buf); fflush(0);
	fprintf(fres, "%d\n", ndata);
	for(i=0; i<ndata; ++i)
		fprintf(fres, "%lf %lf %lf\n", force[i], AD_sim[i], TD_sim[i]);
	fclose(fres);

	/* compute fitness */
	gsl_vector * diff = gsl_vector_calloc(2*ndata);
	gsl_matrix * cov = gsl_matrix_calloc(2*ndata, 2*ndata);
	gsl_matrix * cov_inv = gsl_matrix_calloc(2*ndata, 2*ndata);
	gsl_permutation * perm = gsl_permutation_alloc(2*ndata);
	gsl_vector * temp = gsl_vector_calloc(2*ndata);
	int signum;
	double log_cov_part;

	/* covariance matrix, AD-AD and TD-TD part */
	for(i=0; i<ndata; ++i)
	{
		gsl_vector_set(diff, i, AD_exp[i]-AD_sim[i]);
		gsl_vector_set(diff, ndata+i, TD_exp[i]-TD_sim[i]);
		gsl_matrix_set(cov, i, i, sigma2);
		gsl_matrix_set(cov, ndata+i, ndata+i, sigma2);
		for(j=i+1; j<ndata; ++j)
		{
			double curr = sigma2*exp(-mult_len*pow(force[i]-force[j],2));
			/* AD-AD */
			gsl_matrix_set(cov, i, j, curr);
			gsl_matrix_set(cov, j, i, curr);
			/* TD-TD */
			gsl_matrix_set(cov, ndata+i, ndata+j, curr);
			gsl_matrix_set(cov, ndata+j, ndata+i, curr);
		}
	}

	/* Note: covariance matrix now has a block-diagonal structure, no correlation between AD and TD points */
	/* Note: can speed-up the calculation of the inverse by inverting only one block */

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

	free(force);
	free(AD_exp);
	free(AD_sim);
	free(TD_exp);
	free(TD_sim);

	return res;
}

