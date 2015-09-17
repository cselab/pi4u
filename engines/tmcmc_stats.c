/*
 *  tmcmc_stats.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 1/1/14.
 *  Copyright 2014 ETH Zurich. All rights reserved.
 *
 */

#include <stdio.h>
#include <math.h>
#include "engine_tmcmc.h"
#include <time.h>

int display = 0;

/** OBJLOGP FUNCTION **/
double Objlogp(double x, double *fj, int fn, double pj, double tol)
{
	int i;
	double fjmax = compute_max(fj, fn);

	double weight[fn];
	for (i = 0; i < fn; i++)
		weight[i] = exp((fj[i]-fjmax)*(x-pj));

	double sum_weight = compute_sum(weight, fn);

	double q[fn];
	for (i = 0; i < fn; i++)
		q[i] = weight[i]/sum_weight;

	double mean_q = compute_mean(q, fn);
	double std_q = compute_std(q, fn, mean_q);

	double CoefVar = pow(std_q/mean_q-tol, 2);	/* result */

	return CoefVar;
}

typedef struct fparam_s {
	double *fj;
	int     fn;
	double  pj;
	double  tol;
} fparam_t;

fparam_t *sfp;

double Objlogp_s(double *x, int n)
{
	double *fj = sfp->fj;
	int fn = sfp->fn;
	double pj = sfp->pj;
	double tol = sfp->tol;

	return Objlogp(x[0], fj, fn, pj, tol);
}

double Objlogp_gsl(double x, void *param)
{
	fparam_t *fp = (fparam_t *) param;

	double *fj = fp->fj;
	int fn = fp->fn;
	double pj = fp->pj;
	double tol = fp->tol;

	double res = Objlogp(x, fj, fn, pj, tol);
/*	printf("Objlogp(%lf)=%lf\n", x, res); */
	return res;
}

double Objlogp_gsl2(const gsl_vector *v, void *param)
{
	double x;
	x = gsl_vector_get(v, 0);

	return Objlogp_gsl(x, param);
}

/*** OPTIMIZATION ***/
int fminsearch(double *fj, int fn, double pj, double tol, double *xmin, double *fmin)
{
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minex_func;
	int conv = 0;

	size_t iter = 0, max_iter = data.options.MaxIter;	/* USER input*/
	double Tol = data.options.Tol;
	int Display = data.options.Display;
	double Step = data.options.Step;
	int status;
	double size;

	fparam_t fp;

	fp.fj = fj; fp.fn = fn; fp.pj = pj; fp.tol = tol;

	/* Starting point */
	x = gsl_vector_alloc (1);
	gsl_vector_set (x, 0, pj);

	/* Set initial step sizes to 1 */
	ss = gsl_vector_alloc (1);
	gsl_vector_set_all (ss, Step); /* input?*/

       /* Initialize method and iterate */
	minex_func.n = 1;
	minex_func.f = Objlogp_gsl2;
	minex_func.params = &fp;

       s = gsl_multimin_fminimizer_alloc (T, 1);
       gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, Tol);

		if (status == GSL_SUCCESS) {
			conv = 1;
			if (Display)
				printf ("converged to minimum at\n");
		}

		if (Display)
			printf ("%3ld %.16lf f() = %.16f size = %.16f\n",
				iter, gsl_vector_get (s->x, 0), s->fval, size);

	} while (status == GSL_CONTINUE && iter < max_iter);

	if (conv) {
		*fmin = s->fval;
		*xmin = gsl_vector_get(s->x, 0);
	} else {
		*fmin = 0;
		*xmin = 0.0;
	}

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return conv;
}

int fmincon(double *fj, int fn, double pj, double tol, double *xmin, double *fmin)
{
	int status;
	int iter = 0, max_iter = data.options.MaxIter;	/* USER input*/
	double Tol = data.options.Tol;
	int Display = data.options.Display;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	double m = 0.5;     /* input */
 	double a = 0.0, b = 2.0;    /* input */
	gsl_function F;
	int conv = 0;
	gsl_vector *x;
	int i;

	x = gsl_vector_alloc (1);

	fparam_t fp;

	fp.fj = fj; fp.fn = fn; fp.pj = pj; fp.tol = tol;

	F.function = Objlogp_gsl;
	F.params = &fp;

	T = gsl_min_fminimizer_goldensection; /*brent;*/
	s = gsl_min_fminimizer_alloc (T);

/*	printf("f(a)=%lf\n", Objlogp_gsl(a, &fp));*/
/*	printf("f(b)=%lf\n", Objlogp_gsl(b, &fp));*/
/*	printf("f(m)=%lf\n", Objlogp_gsl(m, &fp));*/
	double fa = Objlogp_gsl(a, &fp);
	double fb = Objlogp_gsl(b, &fp);
	for (i = 0; i < max_iter; i++) {
		m = a + i*(b-a)/max_iter;
		double fm = Objlogp_gsl(m, &fp);
		if ((fm < fa) && (fm < fb)) break;
	}

	if (i == max_iter) {
		if (Display)
			printf("failed to initialize fmincon!\n");
		return 0;
	}
	else {
		if (Display)
			printf("inited with %d tries\n", i);
	}

	gsl_min_fminimizer_set (s, &F, m, a, b);
/*	printf("bbb\n");*/

	if (Display) {
		printf ("using %s method\n", gsl_min_fminimizer_name (s));
		printf ("%5s [%18s, %18s] %18s %18s\n", "iter", "lower", "upper", "min", "err(est)");
		printf ("%5d [%.16f, %.16f] %.16f %.16f\n", iter, a, b, m, b - a);
	}

	do {
		iter++;
		status = gsl_min_fminimizer_iterate (s);

		m = gsl_min_fminimizer_x_minimum (s);
		a = gsl_min_fminimizer_x_lower (s);
		b = gsl_min_fminimizer_x_upper (s);

		status = gsl_min_test_interval (a, b, Tol, 0.0);
		if (status == GSL_SUCCESS) {
			if (Display)
				printf ("Converged:\n");
			conv = 1;
		}

		if (Display)
			printf ("%5d [%.16f, %.16f]  %.16f %.16f\n",
				iter, a, b, m, b - a);

	} while (status == GSL_CONTINUE && iter < max_iter);

	if (conv) {
		gsl_vector_set (x, 0, m);
		*fmin = Objlogp_gsl(m, &fp);
		*xmin = m;
	} else {
		*fmin = 0;
		*xmin = 0.0;
	}

	gsl_vector_free(x);
	gsl_min_fminimizer_free (s);

	return conv;

}


/*** STATISTICS ***/
#include "posdef.c" 	// peh

void calculate_statistics(double flc[], int n, int nselections, int gen, unsigned int sel[])
{
	/*double pflag = 0;*/
	double tolCOV = data.TolCOV;
	double *CoefVar = runinfo.CoefVar;
	double *p = runinfo.p;
	int *Num = data.Num;
/*	int *currentuniques = runinfo.currentuniques; */
	double *logselection = runinfo.logselection;

	double fmin = 0, xmin = 0;
	int conv = 0;
/*	conv = fmincon(flc, n, p[gen], tolCOV, &xmin, &fmin);*/
	if (display)
		printf("fmincon: conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
	if (!conv) {
		conv = fminsearch(flc, n, p[gen], tolCOV, &xmin, &fmin);
		if (display)
			printf("fminsearch: conv=%d xmin=%.16lf fmin=%.16lf\n", conv, xmin, fmin);
	}

	/* gen: next generation number */
	int j = gen+1;

	p[j] = xmin;
	CoefVar[j] = fmin;

	if (p[j]>1) {
		/*pflag=p[j-1];*/
		p[j] = 1;
		Num[j]=data.LastNum;
	}

	/*print_matrix("p", p, j);*/

/*
	if (p[j]>10) { // data.pstrict
		 data.Rsq = data.Rsqstrict;
	}
*/
	/* Compute weights and normalize*/

	int i;
#if 1
	double flcp[n];
	for (i= 0; i<n; i++)
		flcp[i] = flc[i]*(p[j]-p[j-1]);


	double fjmax= compute_max (flcp,n );
	double weight[n];
	/*PA weight[i] = exp((flc[i]-fjmax)*(p[j]-p[j-1])); 23/06 */
	for (i = 0; i < n; i++)
		weight[i] = exp( flcp[i] - fjmax );

	if (display)
		print_matrix("weight", weight, n);

	double sum_weight = compute_sum(weight, n);

	double q[n];
	for (i = 0; i < n; i++)
		q[i] = weight[i]/sum_weight;

	if (display)
		print_matrix("runinfo_q", q, n);

	/*double sum_q = compute_sum(q, n);*/

	/*logselection[gen] = log(sum_weight/currentuniques[gen])+fjmax*(p[gen+1]-p[gen]); PA definition change for all types of resampling 23/06/15*/
	logselection[gen]= log(sum_weight) + fjmax -log(n);
#else
	double fjmax = compute_max(flc, n);

	double weight[n];
	for (i = 0; i < n; i++)
		weight[i] = exp((flc[i]-fjmax)*(p[j]-p[j-1]));

	if (display)
	print_matrix("weight", weight, n);

	double sum_weight = compute_sum(weight, n);

	double q[n];
	for (i = 0; i < n; i++)
		q[i] = weight[i]/sum_weight;

	if (display)
		print_matrix("runinfo_q", q, n);

	/*double sum_q = compute_sum(q, n);*/

	logselection[gen] = log(sum_weight/currentuniques[gen])+fjmax*(p[gen+1]-p[gen]);
#endif

	if (display)
		print_matrix("logselection", logselection, gen+1);

	double mean_q = compute_mean(q, n);
	double std_q = compute_std(q, n, mean_q);

	CoefVar[gen] = std_q/mean_q;

	if (display)
		print_matrix("CoefVar", CoefVar, gen+1);

	size_t K = n;
	unsigned int N = 1;

	unsigned int samples = n; /*1000;*/
	unsigned int nn[samples];

	for (i = 0; i < samples; i++) sel[i] = 0;

	int k;

	if (nselections == 0) nselections = samples; /* n;*/
	for (k = 0; k < nselections; k++) {

		/*gsl_ran_multinomial (r, K, N, q, nn);*/
		multinomialrand (K, N, q, nn);
		for (i = 0; i < K; i++) sel[i]+=nn[i];
	}

	if (display) {
		printf("\n s = [");
		for (i = 0; i < K; i++) printf("%d ", sel[i]);
		printf("]\n");
	}

	/* compute SS */
	int PROBDIM = data.Nth;

	double mean_of_theta[PROBDIM];

	for (i = 0; i < PROBDIM; i++) {
		mean_of_theta[i] = 0;
		for (j = 0; j < n; j++) mean_of_theta[i]+=curgen_db.entry[j].point[i]*q[j];

		runinfo.meantheta[gen][i] = mean_of_theta[i];
	}

	if (display)
		print_matrix("mean_of_theta", mean_of_theta, PROBDIM);

	double meanv[PROBDIM];
	for (i = 0; i < PROBDIM; i++) {
		meanv[i] = mean_of_theta[i];
	}

	for (i = 0; i < PROBDIM; i++) {
		for (j = 0; j < PROBDIM; j++) {
			double s;
			int k;
			s = 0;
			for (k = 0; k < n; k++) {
				s += q[k]*(curgen_db.entry[k].point[i]-meanv[i])*(curgen_db.entry[k].point[j]-meanv[j]);
			}
			runinfo.SS[i][j] = runinfo.SS[j][i] = s;
		}
	}

#if 1	/* peh:check this */
	{
	int fixed = make_posdef(runinfo.SS[0], PROBDIM, 2);
	if (fixed) {
		printf("WARNING: runinfo.SS was forced to become positive definite\n");
	}
	}
#endif

	if (display)
		print_matrix_2d("runinfo.SS", runinfo.SS, PROBDIM, PROBDIM);
}

double priorpdf(double *theta, int n)
{
	/* peh:check this */
	double res;

	if (data.prior_type == 0)
	{
		/* log-uniform */
		res = 0;
		/*res = 1;*/

		int i;
		for (i = 0; i < n; i++) {
			/*res += log(gsl_ran_flat_pdf(theta[i], data.lowerbound[i], data.upperbound[i])); PA EDIT
			res *= gsl_ran_flat_pdf(theta[i], data.lowerbound[i], data.upperbound[i]); */
			res += -log( data.upperbound[i]- data.lowerbound[i]);
		}
		if (res == 0) return 0;
		/*return log(res); PA EDIT*/
		return res;
	}
	else
	{
		/* gaussian */
		res = logmvnpdf(n, theta, data.prior_mu, data.prior_sigma);
	}
	return res;
}

double posterior(double *theta, int n, double LH)
{
	double res;
/*
	double Prior = priorpdf(theta, n);
	printf("Prior = %lf\n", Prior);

	if (Prior > 0)
		res = LH + log(Prior);	// xxx
	else
		res = LH;
*/

	res = LH;	/* Algorithm fix by PanosA */

	return res;
}
