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
#include <sys/time.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ftw.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <stdlib.h>
#include <torc.h>
#include "gsl_headers.h"
#include "spawner.h"

#define DBG 0

int dbg_display = 1;

const char * from_dir = "/users/lina/UQ/pi4u.git/engines/fitfun_code/mpi-dpd";
const char *data_filename = "diam_vs_force_data_u.txt"; /* _u -- with uncertainty */
const int ndata = 7;

void sighandler(int signo)
{
	printf("spawner: SIGNAL %d was received!\n", signo); fflush(0);
}

void copy_files_to_fitfundir(const char * fitfun_dir)
{
	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	chdir(fitfun_dir);
	copy_file(from_dir, data_filename);
	copy_file(from_dir, "diam_plot.gp");

	chdir(from_dir);
	chdir("../cuda-rbc");
	char cudadir[1024];
	getcwd(cudadir, sizeof(cudadir));

	chdir(fitfun_dir);
	mkdir("cuda-rbc", S_IRWXU);
	chdir("cuda-rbc");
	copy_file(cudadir, "rbc.dat");

	chdir(workdir);
}

void rm_files_from_fitfundir(const char * fitfun_dir)
{
	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	chdir(fitfun_dir);
	rmrf("cuda-rbc");
	rmrf("diam_plot.gp");

	chdir(workdir);
}

void copy_files_to_fitfuntaskdir(const char * fitfuntask_dir)
{
	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	chdir(fitfuntask_dir);
	copy_file(from_dir, "test");
	copy_file(from_dir, "rbcs-ic.txt");
	copy_file(from_dir, "compute_diameters_standalone.py");
	copy_file(from_dir, "plyfile.py");
	copy_file(from_dir, "load_modules.sh");

	chdir(workdir);
}

void rm_files_from_fitfuntaskdir(const char * fitfuntask_dir)
{
	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	chdir(fitfuntask_dir);

	copy_file("ply", "rbcs-0000.ply");
	copy_file("ply", "rbcs-0001.ply");
	copy_file("ply", "rbcs-0099.ply");

	rmrf("ply");
	rmrf("xyz");
	rmrf("test");
	rmrf("rbcs-ic.txt");
	rmrf("compute_diameters_standalone.py");
	rmrf("plyfile.py");
	rmrf("plyfile.pyc");
	rmrf("load_modules.sh");
	rmrf("diag.txt");
	rmrf("output.txt");
	rmrf("diam.txt");

	chdir(workdir);
}

void fitfuntask(double in_tparam[], double *out_tparam,	char fitfun_dir[])
{
	int me = getpid();	/* spawner_id : worker_id */
	int rf;

	double x0		= in_tparam[0];
	double p_F		= in_tparam[1];
	double force_P	= in_tparam[2];

	char line[1024];
	char args[1024];
	char *largv[64];

	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	char fitfuntask_dir[1024];
	sprintf(fitfuntask_dir, "%s/output_force_%.16lf", fitfun_dir, force_P);

	mkdir(fitfuntask_dir, S_IRWXU);
	copy_files_to_fitfuntaskdir(fitfuntask_dir);
	chdir(fitfuntask_dir);	/* enter the new directory */

	rf = fork();
	if (rf < 0)
	{
		printf("spawner(%d): fork failed!!!!\n", me);
		fflush(0);
	}
	if (rf == 0)
	{
		const double Fscale = 4.2e-14; // units conversion from master
		double force_F = force_P / Fscale;

		strcpy(line, "");
		strcat(line, "./test 1 1 1 -rbcs -tend=100 -steps_per_dump=1000");

		strcpy(args, "");
		sprintf(args, " -stretchingForce=%lf -RBCtotArea=124 -RBCtotVolume=90 -RBCka=4900 -RBCkb=40 -RBCkd=100 -RBCkv=5000 -RBCgammaC=30 -RBCx0=%lf -RBCp=%f", force_F, x0, p_F);
		strcat(line, args);
#if DBG
		printf("Running command for diam vs force %s\n", line);
#endif

		char out_file[1024];
		sprintf(out_file, "%s/output.txt", fitfuntask_dir);

		int fd = open(out_file, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient

		parse(line, largv);	/* prepare argv */
		execvp(*largv, largv);
	}
	waitpid(rf, NULL, 0);
	sync();

	char out_file_post[1024];
	sprintf(out_file_post, "%s/diam.txt", fitfuntask_dir);

	/* get the AD, TD */
	rf = fork();
	if (rf < 0)
	{
		printf("spawner(%d): fork failed!!!!\n", me);
		fflush(0);
	}
	if (rf == 0)
	{
		strcpy(line, "");
		strcat(line, "python compute_diameters_standalone.py --ply_filedir=ply");

		int fd = open(out_file_post, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
		dup2(fd, 1);	// make stdout go to file
		dup2(fd, 2);	// make stderr go to file
		close(fd);		// fd no longer needed - the dup'ed handles are sufficient

		parse(line, largv);	/* prepare argv */
		execvp(*largv, largv);
	}
	waitpid(rf, NULL, 0);
	sync();

	FILE * fp = fopen(out_file_post, "r");
	if(!fp)
		printf("Cannot open file %s\n", out_file_post);
	else
	{
		fscanf(fp, "%lf %lf %lf %lf", &out_tparam[0], &out_tparam[1], &out_tparam[2], &out_tparam[3]);
#if DBG
		printf("%lf %lf %lf %lf %lf\n", force_P, out_tparam[0], out_tparam[1], out_tparam[2], out_tparam[3]);
#endif

		fclose(fp);
	}

	chdir(workdir);

#if !DBG
	rm_files_from_fitfuntaskdir(fitfuntask_dir);
#endif
}

double fitfun(double *input, int n, void *output, int *info)
{
	const char *res_filename = "diam_vs_force_result.txt";
	double res = -1e6; /* return in case of failure */
	int me = getpid();	/* spawner_id : worker_id */
	int i;

	double x0		= input[0];
	double p_F		= input[1] * 1e-3; // descale
	double sigmaAD	= input[2];
	double len		= input[3] * 1e2; // descale

	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	char fitfun_dir[256];
	sprintf(fitfun_dir, "%s/fitfun.%d.%d.%d.%d", workdir, info[0], info[1], info[2], info[3]);
	mkdir(fitfun_dir, S_IRWXU);

	if (dbg_display)
	{
		printf("spawner(%d): running task %s with params (", me, fitfun_dir);

		for (i = 0; i < n-1; i++)
			printf("%.16lf, ", input[i]);
		printf("%.16lf)\n", input[i]);

		fflush(0);
	}

	double t0, t1;
	t0 = my_gettime();
	/* 1. PREPROCESSING PHASE - APPLICATION SPECIFIC */
	/* 1a. Write parameters to the file */
	char parfile[1024];
	sprintf(parfile, "%s/params.dat", fitfun_dir);
	FILE *finp = fopen(parfile, "w");
	for (i = 0; i < n; i++)
		fprintf(finp, "%.16lf\n", input[i]);
	fclose(finp);

	/* 1b. Copy some application specific input files*/
	copy_files_to_fitfundir(fitfun_dir);

	/* 2. EXECUTION OF THE FUNCTION EVALUATION (SIMULATION) SOFTWARE */
	/* 2a. Read data from file */
	char datafile[1024];
	sprintf(datafile, "%s/%s", fitfun_dir, data_filename);
	FILE * fdata = fopen(datafile, "r");
	if(!fdata)
	{
		printf("Cannot open file %s\n", datafile);
		fflush(0);
		return res;
	}

	/* allocate arrays */
	double * force_P		= malloc(ndata*sizeof(double));
	double * ADexp			= malloc(ndata*sizeof(double));
	double * TDexp			= malloc(ndata*sizeof(double));
	double * ADexp_u		= malloc(ndata*sizeof(double));
	double * TDexp_u		= malloc(ndata*sizeof(double));
	gsl_vector * diff		= gsl_vector_calloc(ndata);
	gsl_matrix * cov		= gsl_matrix_calloc(ndata, ndata);
	gsl_matrix * cov_inv	= gsl_matrix_calloc(ndata, ndata);
	gsl_permutation * perm	= gsl_permutation_alloc(ndata);
	gsl_vector * temp		= gsl_vector_calloc(ndata);

	for(i=0; i<ndata; ++i)
		fscanf(fdata, "%lf %lf %lf %lf %lf\n", &force_P[i], &ADexp[i], &ADexp_u[i],
			&TDexp[i], &TDexp_u[i]);
	fclose(fdata);

	/* 2b. Compute stretching vs force */
	double * out_tparam = malloc(4*ndata*sizeof(double)); /* ADsim[i] ADsim_u[i] TDsim[i] TDsim_u[i] */
	double in_tparam[3];
	for(i=0; i<ndata; ++i)
	{
		/* Run the simulation with current force */
		in_tparam[0] = x0;
		in_tparam[1] = p_F;
		in_tparam[2] = force_P[i] * 1e-12; // descale
		torc_create(-1, fitfuntask, 3,
			3,					MPI_DOUBLE, CALL_BY_COP,
			4,					MPI_DOUBLE,	CALL_BY_RES,
			strlen(fitfun_dir),	MPI_INT,	CALL_BY_COP, /* TODO: check with Panos if MPI_CHAR is supported*/
			in_tparam, out_tparam+4*i, fitfun_dir);
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

	/* 3. POSTPROCESSING PHASE - APPLICATION SPECIFIC */
	/* write results to file */
	char resfile[1024];
	sprintf(resfile, "%s/%s", fitfun_dir, res_filename);
	FILE * fres = fopen(resfile, "w");
	if(!fres)
	{
		printf("Cannot open file %s\n", resfile);
		fflush(0);
	}
	for(i=0; i<ndata; ++i)
		fprintf(fres, "%lf %lf %lf %lf %lf\n", force_P[i], out_tparam[4*i], out_tparam[4*i+1],
			out_tparam[4*i+2], out_tparam[4*i+3]);
	fclose(fres);

	/* plot fit */
	int rf = fork();
	if (rf < 0)
	{
		printf("spawner(%d): fork failed!!!!\n", me);
		fflush(0);
	}
	if (rf == 0)
	{
		char line[1024];
		char *largv[64];

		chdir(fitfun_dir);

		strcpy(line, "");
		strcat(line, "gnuplot diam_plot.gp");

		parse(line, largv);	/* prepare argv */
		execvp(*largv, largv);
	}
	waitpid(rf, NULL, 0);
	sync();

	chdir(fitfun_dir);
	char plotname[1024];
	sprintf(plotname, "diam_vs_force_%.2e_%.2e.png", x0, p_F);
	rename("diam_vs_force.png", plotname);
	chdir(workdir);

	/* check for NaN and Inf and zero */
	for(i=0; i<4*ndata; ++i)
		if (isnan(out_tparam[i]) || isinf(out_tparam[i]) || out_tparam[i] < 1e-6)
		{
			res = -1e12;
			goto finalize;
		}

	int signum;
	double log_cov_part, s;

	double sigmaAD2		= sigmaAD*sigmaAD;
	double mult_len		= 0.5/(len*len);

	/* covariance matrix, AD-AD and TD-TD part */
	for(i=0; i<ndata; ++i)
	{
		double ADexp_i = ADexp[i];
		double ADsim_i = out_tparam[4*i];

		/* theta-mu */
		gsl_vector_set(diff, i, ADexp_i-ADsim_i);

		double ADexp_u_i = ADexp_u[i];
		double ADsim_u_i = out_tparam[4*i+1];

		/* AD-AD diag */
		gsl_matrix_set(cov, i, i, sigmaAD2 + ADexp_u_i*ADexp_u_i + ADsim_u_i*ADsim_u_i);

		// model error part
		int j;
		for(j=i+1; j<ndata; ++j)
		{
			double fij = abs(force_P[i]-force_P[j]);
			if (fij < len)
			{
				double mult = exp(-mult_len*fij*fij);

				/* AD-AD */
				s = mult*sigmaAD2;
				gsl_matrix_set(cov, i, j, s);
				gsl_matrix_set(cov, j, i, s);
			}
		}
	}

	/* Note: if the covariance matrix has a block-diagonal structure, one can speed-up the calculation of the inverse by inverting block by block */

	gsl_set_error_handler_off();

	int status;
	status = gsl_linalg_LU_decomp(cov, perm, &signum); /* now cov contains its LU decomposition */
	if(status != GSL_SUCCESS)
	{
		fprintf (stderr, "fitfun.c: LU decomposition failed\n");
		fflush(0);
		goto finalize;
	}

	status = gsl_linalg_LU_invert(cov, perm, cov_inv); /* cov_inv contains inverse of cov */
	if(status != GSL_SUCCESS)
	{
		fprintf (stderr, "fitfun.c: LU invertion failed\n");
		fflush(0);
		goto finalize;
	}

	double log_det_cov = gsl_linalg_LU_lndet(cov);
	gsl_blas_dgemv(CblasNoTrans, 1.0, cov_inv, diff, 0.0, temp);
	gsl_blas_ddot(diff, temp, &log_cov_part);

	res = -0.5 * ( ndata*log(2*M_PI) + log_det_cov + log_cov_part );

	if (isnan(res) || isinf(res))
		res = -1e12;

finalize:

	if (dbg_display)
	{
		printf("task(%d.%d.%d.%d): ", info[0], info[1], info[2], info[3]);

		for (i=0; i<n-1; i++)
			printf("%.16lf, ", input[i]);
		printf("%.16lf ", input[i]);

		printf("= %lf ", res);

		printf("in %lf secs\n", t1-t0);

		fflush(0);
	}

#if !DBG
	rm_files_from_fitfundir(fitfun_dir);
#endif

	/* free memory */
	gsl_vector_free(diff);
	gsl_matrix_free(cov);
	gsl_matrix_free(cov_inv);
	gsl_permutation_free(perm);
	gsl_vector_free(temp);

	free(force_P);
	free(ADexp);
	free(ADexp_u);
	free(TDexp);
	free(TDexp_u);
	free(out_tparam);

	return res;
}

