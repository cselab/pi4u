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

#define FPF (1e-12/7.2237e-14)
#define DBG 0
//#define CONTACTFORCES

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
	int me = getpid();	/* spanwer_id : worker_id */
	int rf;

	double lmax		= in_tparam[0];
	double p		= in_tparam[1];
	double force	= in_tparam[2]*FPF;

	char line[1024];
	char args[1024];
	char *largv[64];

	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	char fitfuntask_dir[1024];
	sprintf(fitfuntask_dir, "%s/output_force_%.16lf", fitfun_dir, force);

	mkdir(fitfuntask_dir, S_IRWXU);
	copy_files_to_fitfuntaskdir(fitfuntask_dir);
	chdir(fitfuntask_dir);	/* enter to the new directory */

	rf = fork();
	if (rf < 0)
	{
		printf("spanwer(%d): fork failed!!!!\n", me);
		fflush(0);
	}
	if (rf == 0)
	{
		double l0 = 0.559432;
		double x0 = l0/lmax;
		int q = 1;
		double A0 = sqrt(3)*l0*l0*0.25;
		double kBT = 0.02883048889;
		double cq = sqrt(3)*pow(A0,q+1)*kBT*(4*x0*x0 - 9*x0 + 6) / (4*p*q*lmax*(1-x0)*(1-x0));

		strcpy(line, "");
#ifdef CONTACTFORCES
		strcat(line, "./test 1 1 1 -rbcs -tend=100 -steps_per_dump=1000 -contactforces");
#else
		strcat(line, "./test 1 1 1 -rbcs -tend=100 -steps_per_dump=1000");
#endif
		strcpy(args, "");
		sprintf(args, " -RBClmax=%lf -RBCp=%lf -RBCcq=%lf -stretchingForce=%lf", lmax, p, cq, force);
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
		printf("spanwer(%d): fork failed!!!!\n", me);
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
		printf("%lf %lf %lf %lf %lf\n", force, out_tparam[0], out_tparam[1], out_tparam[2], out_tparam[3]);
#endif

		fclose(fp);
	}

	chdir(workdir);

	rm_files_from_fitfuntaskdir(fitfuntask_dir);
}

double fitfun(double *input, int n, void *output, int *info)
{
	const char *res_filename = "diam_vs_force_result.txt";
	double res = -1e6; /* return in case of failure */
	int me = getpid();	/* spanwer_id : worker_id */
	int i, j;

	double lmax		= input[0];
	double p		= input[1];
	double sigmaAD	= input[2];
	double sigmaTD	= input[3];

	char workdir[1024];
	getcwd(workdir, sizeof(workdir));

	char fitfun_dir[256];
	sprintf(fitfun_dir, "%s/fitfun.%d.%d.%d.%d", workdir, info[0], info[1], info[2], info[3]);
	mkdir(fitfun_dir, S_IRWXU);

	if (dbg_display)
	{
		printf("spanwer(%d): running task %s with params (", me, fitfun_dir);

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
	double * force = malloc(ndata*sizeof(double));
	double * ADexp = malloc(ndata*sizeof(double));
	double * TDexp = malloc(ndata*sizeof(double));
	double * ADexp_u = malloc(ndata*sizeof(double));
	double * TDexp_u = malloc(ndata*sizeof(double));
	for(i=0; i<ndata; ++i)
		fscanf(fdata, "%lf %lf %lf %lf %lf\n", &force[i], &ADexp[i], &ADexp_u[i],
			&TDexp[i], &TDexp_u[i]);
	fclose(fdata);

	/* 2b. Compute stretching vs force */
	double * out_tparam = malloc(4*ndata*sizeof(double)); /* ADsim[i] ADsim_u[i] TDsim[i] TDsim_u[i] */
	double in_tparam[3];
	for(i=0; i<ndata; ++i)
	{
		/* Run the simulation with current force */
		in_tparam[0] = lmax;
		in_tparam[1] = p;
		in_tparam[2] = force[i];
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
		
	if (dbg_display)
		printf("spanwer(%d): simcode_time=%lf secs\n", me, t1-t0);fflush(0);

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
		fprintf(fres, "%lf %lf %lf %lf %lf\n", force[i], out_tparam[4*i], out_tparam[4*i+1],
			out_tparam[4*i+2], out_tparam[4*i+3]);
	fclose(fres);

	/* plot fit */
	int rf = fork();
	if (rf < 0)
	{
		printf("spanwer(%d): fork failed!!!!\n", me);
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
	sprintf(plotname, "diam_vs_force_%.2e_%.2e.png", lmax, p);
	rename("diam_vs_force.png", plotname);
	chdir(workdir);

	/* compute fitness */
	gsl_vector * diff		= gsl_vector_calloc(2*ndata);
	gsl_matrix * cov		= gsl_matrix_calloc(2*ndata, 2*ndata);
	gsl_matrix * cov_inv	= gsl_matrix_calloc(2*ndata, 2*ndata);
	gsl_permutation * perm	= gsl_permutation_alloc(2*ndata);
	gsl_vector * temp		= gsl_vector_calloc(2*ndata);

	int signum;
	double log_cov_part, s;

	double sigmaAD2		= sigmaAD*sigmaAD;
	double sigmaTD2		= sigmaTD*sigmaTD;
	double sigmaADTD	= sigmaAD*sigmaTD;
	double len			= 30;
	double mult_len		= 0.5/(len*len);

	/* covariance matrix, AD-AD and TD-TD part */
	for(i=0; i<ndata; ++i)
	{
		int ADi = i;
		int TDi = i + ndata;

		double ADexp_i = ADexp[i];
		double TDexp_i = TDexp[i];
		double ADsim_i = out_tparam[4*i];
		double TDsim_i = out_tparam[4*i+2];

		/* theta-mu */
		gsl_vector_set(diff, ADi, ADexp_i-ADsim_i);
		gsl_vector_set(diff, TDi, TDexp_i-TDsim_i);

		double ADexp_u_i = ADexp_u[i];
		double TDexp_u_i = TDexp_u[i];
		double ADsim_u_i = out_tparam[4*i+1];
		double TDsim_u_i = out_tparam[4*i+3];

		/* AD-AD diag */
		gsl_matrix_set(cov, ADi, ADi, sigmaAD2 + ADexp_u_i*ADexp_u_i + ADsim_u_i*ADsim_u_i);

		/* TD-TD diag */
		gsl_matrix_set(cov, TDi, TDi, sigmaTD2 + TDexp_u_i*TDexp_u_i + TDsim_u_i*TDsim_u_i);

		/* AD-TD diag */
		s = sigmaADTD + ADexp_u[i]*TDexp_u[i] + ADsim_u_i*TDsim_u_i;
		gsl_matrix_set(cov, i, TDi, s);
		gsl_matrix_set(cov, TDi, i, s);

		for(j=i+1; j<ndata; ++j)
		{
			int ADj = j;
			int TDj = j + ndata;
	
			double ADexp_u_j = ADexp_u[j];
			double TDexp_u_j = TDexp_u[j];

			double ADsim_u_j = out_tparam[4*j+1];
			double TDsim_u_j = out_tparam[4*j+3];

			double fij = abs(force[i]-force[j]);
			if (fij < len)
			{
				double mult = exp(-mult_len*fij*fij);
	
				/* AD-AD */
				s = mult * (sigmaAD2 + ADexp_u_i*ADexp_u_j + ADsim_u_i*ADsim_u_j);
				gsl_matrix_set(cov, ADi, ADj, s);
				gsl_matrix_set(cov, ADj, ADi, s);
	
				/* TD-TD */
				s = mult * (sigmaTD2 + TDexp_u_i*TDexp_u_j + TDsim_u_i*TDsim_u_j);
				gsl_matrix_set(cov, TDi, TDj, s);
				gsl_matrix_set(cov, TDj, TDi, s);
			}
		}
	}

	/* Note: if the covariance matrix has a block-diagonal structure, one can speed-up the calculation of the inverse by inverting block by block */

	gsl_set_error_handler_off();

	int status;
	status = gsl_linalg_LU_decomp(cov, perm, &signum); /* now cov contains its LU decomposition */
	if(status != GSL_SUCCESS)
	{
		fprintf (stderr, "LU failed, gsl_errno=%d\n", status);
		fflush(0);
		goto done;
	}

	status = gsl_linalg_LU_invert(cov, perm, cov_inv); /* cov_inv contains inverse of cov */
	if(status != GSL_SUCCESS)
	{
		fprintf (stderr, "Invertion failed, gsl_errno=%d\n", status);
		fflush(0);
		goto done;
	}

	double log_det_cov = gsl_linalg_LU_lndet(cov);
	gsl_blas_dgemv(CblasNoTrans, 1.0, cov_inv, diff, 0.0, temp);
	gsl_blas_ddot(diff, temp, &log_cov_part);

	res = -0.5 * ( ndata*log(2*M_PI) + log_det_cov + log_cov_part );

done:

	if(dbg_display)
		printf("log-likelihood: %.16f\n", res);

	rm_files_from_fitfundir(fitfun_dir);

	gsl_vector_free(diff);
	gsl_matrix_free(cov);
	gsl_matrix_free(cov_inv);
	gsl_permutation_free(perm);
	gsl_vector_free(temp);

	free(force);
	free(ADexp);
	free(ADexp_u);
	free(TDexp);
	free(TDexp_u);
	free(out_tparam);

	return res;
}

