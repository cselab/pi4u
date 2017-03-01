#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
// testing MCMC with Gaussian target
// high dimensional Gaussian target with positivity constraint

int compare (const void * a, const void * b)
{
	return (*(double*)a - *(double*)b);
}

double sign(double a)
{
	if (a > 0) return +1;
	if (a < 0) return -1;
	return 0;
}

double norm(double *a, int n)
{
	double s = 0;
	for (int i = 0; i < n; i++)
		s += pow(a[i],2.0);

	return sqrt(s);
}

void covcond(double c, double *a, double *Sig, double *Lam, int n)
{
	//%COVCOND covariance matrix with given condition number
	//% [Sig,Lam] = covcond(C,A) generates covariance matrix and its
	//% inverse with given cond number C and first direction A.

	//% create orthogonal basis z, with 1 direction given by 'a'

//	e     = sort(1./linspace(c,1,length(a)));
	double e[n];

	for (int i = 0; i < n; i++)
	{
		e[i] = c + i*((1-c)/(n-1));
		e[i] = 1.0/e[i];
#if DEBUG
		printf("e[%d] = %f\n", i, e[i]);
#endif
	}

	qsort (e, n, sizeof(double), compare);
#if DEBUG
	for (int i = 0; i < n; i++)
	{
		printf("e[%d] = %f\n", i, e[i]);		
	}

	for (int i = 0; i < n; i++)
	{
		printf("a[%d] = %f\n", i, a[i]);		
	}
#endif

//	a(1)  = a(1) + sign(a(1)) * norm(a);  %the Householder trick
	a[0] = a[0] + sign(a[0]) * norm(a, n);

	for (int i = 0; i < n; i++)
	{
		printf("a[%d] = %f\n", i, a[i]);		
	}


//	z     = eye(length(a)) - 2.0/norm(a)^2*a*a';
	double norma2 = pow(norm(a, n), 2);
	double z[n][n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				z[i][j] = 1.0 + (-2.0/norma2)*a[i]*a[j];
			else
				z[i][j] = (-2.0/norma2)*a[i]*a[j];

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("z[%d][%d] = %f\n", i, j, z[i][j]);		
	}
#endif

	double S[n][n], L[n][n];
	double ze[n][n], zt[n][n];

//	Sig   = z * diag(e) * z' ;              % target covariance

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ze[i][j] = z[i][j]*e[j];
			zt[i][j] = z[j][i];
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("ze[%d][%d] = %f\n", i, j, ze[i][j]);		
	}
#endif

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double s = 0.0;
			for (int k = 0; k < n; k++)
				s += ze[i][k]*zt[k][j];
			S[i][j] = s;
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("S[%d][%d] = %f\n", i, j, S[i][j]);		
	}
#endif

//	Lam   = z * inv(diag(e)) * z';          % and its inverse
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ze[i][j] = z[i][j]*(1.0/e[j]);
			zt[i][j] = z[j][i];
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("ze[%d][%d] = %f\n", i, j, ze[i][j]);		
	}
#endif

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double s = 0.0;
			for (int k = 0; k < n; k++)
				s += ze[i][k]*zt[k][j];
			L[i][j] = s;
		}
	}

#if DEBUG
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			printf("L[%d][%d] = %f\n", i, j, L[i][j]);		
	}
#endif

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Sig[i*n+j] = S[i][j];
			Lam[i*n+j] = L[i][j];
		}
	}
}

//1x2 2x2 = 1x2
//
//[a b]	[1  2]
//	[3  4]  
//[x y] [a b]

double *d_mu;
double *d_Lam;

#include "fitfun.c"

double ssfun(double *x, int n)
{
	double s = -fitfun(x, n, NULL, NULL);
#if DEBUG
	printf("s = %f\n", s);
#endif
	return s;
}


double priorfun(double *x, int n)
{
	return 0;
}

struct _params
{
	double *par0;		// = mu+0.1; // initial value
	double sigma2;		// = 1;
	double n0;		// = -1;
	double *lbounds;	// = -inf
	double *ubounds;	// = +inf
	char filename[256];	// = "chain.txt"

	int n;			// = 0;
} params;


void init_params()
{
	params.par0 = NULL;

	params.sigma2 = 1;
	params.n0 = -1;

	params.lbounds = NULL;
	params.ubounds = NULL;
	strcpy(params.filename, "chain.txt");

	params.n = 0;
}

struct _options
{
	int nsimu;		//  =nsimu;
	int adaptint;		// = 500;
	double *qcov;     	//= Sig.*2.4^2./npar;
	double drscale; 	// = drscale;
	double adascale;	// = adascale; // default is 2.4/sqrt(npar) ;

	int printint;		// = 2000;
	double qcovadj;		// = 1e-5;

	int verbosity;		// = 1;
} options;

void init_options()
{
	options.qcovadj = 1e-5;
	options.verbosity = 1;
}


#include "dramrun.c"

int main(int argc, char *argv[])
{
	gsl_rand_init(1);

	//addpath([pwd,filesep,'utils']);

	const int npar     = 2;     // dimension of the target
	double drscale  = 20;     // DR shrink factor
	double adascale = 2.4/sqrt(npar); // scale for adaptation
	int nsimu    = 1e4;    // number of simulations
	int pos = 0;           // positivity?

	// with positivity, you need nsimu to be at least
	// npar ->  nsimu 
	// 10   ->  10 000
	// 20   ->  500 000
	// 100  ->  5 000 000 ?

	// create target covariance matrix, with 1. direction given by 'a'
	// and cond number by 'c'
	double c = 10;  // cond number of the target covariance 

	//a = ones(npar,1); // 1. direction
	double *a = malloc(npar*sizeof(double));
	for (int i = 0; i < npar; i++)
		a[i] = 1;

	double Sig[npar*npar], Lam[npar*npar];
	covcond(c, a, Sig, Lam, npar);
	//[Sig,Lam] = covcond(c,a) // covariance and its inverse

	//mu = zeros(1,npar);       // center point
	double *mu = malloc(npar*sizeof(double));
	for (int i = 0; i < npar; i++)
		mu[i] = 0;

	// create input arguments for the dramrun function
	//clear model params options

	// sum of squares function
	//model.ssfun    = inline('(x-d.mu)*d.Lam*(x-d.mu)''','x','d');

	init_params();
	params.par0    = malloc(npar*sizeof(double));
	for (int i = 0; i < npar; i++)
		params.par0[i] = mu[i]+0.1;	// initial value

	params.lbounds = malloc(npar*sizeof(double));
	for (int i = 0; i < npar; i++)
		params.lbounds[i] = -1e10;	// -inf

	params.ubounds = malloc(npar*sizeof(double));
	for (int i = 0; i < npar; i++)
		params.ubounds[i] = +1e10;	// -inf

	// positivity
	if (pos)
	{
	//	params.bounds = (ones(npar,2)*diag([0,Inf]));
		for (int i = 0; i < npar; i++)
			params.lbounds[i] = 0;
	}

	strcpy(params.filename, "generaltest_chain.txt");

	// arguments for ssfun are in data
	//data  = struct('mu',mu,'Lam',Lam);
	d_mu = mu;
	d_Lam = Lam;

	//double x[2] = {0.5, 0.2};
	//double s = ssfun(x, npar);
	//printf("s = %f\n", s);


	init_options();
	options.nsimu    = nsimu;
	options.adaptint = 500;

	//options.qcov     = Sig.*2.4^2./npar;
	options.qcov = malloc(npar*npar*sizeof(double));
	for (int i = 0; i < npar; i++)
	for (int j = 0; j < npar; j++)
		options.qcov[i*npar+j] = Sig[i*npar+j]*pow(2.4,2)/npar;

	options.drscale  = drscale;
	options.adascale = adascale; // default is 2.4/sqrt(npar) ;
	options.printint = 2000;

	for (int i = 0; i < npar; i++)
		for (int j = 0; j < npar; j++)
			printf("options.qcov[%d][%d]=%f\n", i, j, options.qcov[i*npar+j]);


	dramrun(npar);

//	[results,chain] = dramrun(model,data,params,options);


	return 0;
}
