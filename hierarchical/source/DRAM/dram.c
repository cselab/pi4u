#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fitfun.h"

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
	double e[n];

	for (int i = 0; i < n; i++)
	{
		e[i] = c + i*((1-c)/(n-1));
		e[i] = 1.0/e[i];
	}

	qsort (e, n, sizeof(double), compare);
	a[0] = a[0] + sign(a[0]) * norm(a, n);

	for (int i = 0; i < n; i++)
	{
		printf("a[%d] = %f\n", i, a[i]);		
	}

	double norma2 = pow(norm(a, n), 2);
	double z[n][n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i == j)
				z[i][j] = 1.0 + (-2.0/norma2)*a[i]*a[j];
			else
				z[i][j] = (-2.0/norma2)*a[i]*a[j];

	double S[n][n], L[n][n];
	double ze[n][n], zt[n][n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ze[i][j] = z[i][j]*e[j];
			zt[i][j] = z[j][i];
		}
	}

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


	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			ze[i][j] = z[i][j]*(1.0/e[j]);
			zt[i][j] = z[j][i];
		}
	}


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


	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			Sig[i*n+j] = S[i][j];
			Lam[i*n+j] = L[i][j];
		}
	}
}




double *d_mu;
double *d_Lam;



double ssfun(double *x, int n)
{
	double s = -fitfun(x, n, NULL, NULL);
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


struct _input
{
    	int 	npar;    	/* dimension of the target */
    	double 	drscale;    	/* DR shrink factor */
    	double 	adascale;    	/* scale for adaptation */
    	int 	nsimu;    	/* number of simulations */
    	int 	pos;    	/* positivity [0/1] */
	double	c;		/* cond number of the target covariance */ 
	double 	*lbounds;
	double 	*ubounds;

    	/* with positivity, `nsimu' needs to be at least:
    	   npar ->  nsimu 
    	   10   ->  10 000
    	   20   ->  500 000
    	   100  ->  5 000 000 ? */
} input;



void set_inputs()
{
	/* DEFAULT VALUES */
    	input.npar 	= 4;
    	input.drscale 	= 20;
    	input.adascale 	= 2.4/sqrt(input.npar);
    	input.nsimu 	= 1e5;
    	input.pos 	= 0;
    	input.c 	= 10;
	double 	lb	= -1e10;
	double 	ub	= +1e10;


	/* USER-DEFINED VALUES */
    	FILE *f = fopen("dram.par", "r");
    	if (f == NULL) {
	    printf("No dram.par file available.\n");
    	    exit(1);
    	}

    	char line[256];

    	int line_no = 0;
    	while (fgets(line, 256, f)!= NULL) {
    	    line_no++;
    	    if ((line[0] == '#')||(strlen(line)==0)) {
    	        continue;
    	    }

    	    if (strstr(line, "Npar")) {
    	        sscanf(line, "%*s %d", &input.npar);
    	    }
    	    else if (strstr(line, "drscale")) {
    	        sscanf(line, "%*s %lf", &input.drscale);
    	    }
    	    else if (strstr(line, "adascale")) {
    	        sscanf(line, "%*s %lf", &input.adascale);
    	    }
    	    else if (strstr(line, "nsimu")) {
    	        sscanf(line, "%*s %d", &input.nsimu);
    	    }
    	    else if (strstr(line, "pos")) {
    	        sscanf(line, "%*s %d", &input.pos);
    	    }
    	    else if (strstr(line, "CondNum")) {
    	        sscanf(line, "%*s %lf", &input.c);
    	    }
    	}


	/* Read the lower and upper bounds for each parameter */
    	input.lbounds = (double *)malloc(input.npar*sizeof(double));
    	input.ubounds = (double *)malloc(input.npar*sizeof(double));

	rewind(f);
	line_no = 0;
	int found;

	for (int i = 0; i < input.npar; i++) {
	    found = 0;
	    while (fgets(line, 256, f)!= NULL) {
	        line_no++;
	
	        if ((line[0] == '#')||(strlen(line)==0)) continue;
	
	        char bound[16];
	        sprintf(bound, "B%d", i);
	        if (strstr(line, bound) != NULL) {
	            sscanf(line, "%*s %lf %lf", &input.lbounds[i], &input.ubounds[i]);
	            found = 1;
	            break;
	        }
	    }
	    if (!found) {
	    	input.lbounds[i] = lb;
	    	input.ubounds[i] = ub;
	    }
	}
}



void print_inputs()
{
	/* Print all input (user-defined) parameters */
	printf("\n***************************\n");
	printf("Running with parameters:\n");
	printf(" > npar = %d\n", input.npar);
	printf(" > drscale = %f\n", input.drscale);
	printf(" > adascale = %f\n", input.adascale);
	printf(" > nsimu = %d\n", input.nsimu);
	printf(" > pos = %d\n", input.pos);
	printf(" > c = %f\n", input.c);
	printf(" > Bounds:\n");
	for (int i = 0; i < input.npar; i++)
		printf("  %e %e\n", input.lbounds[i], input.ubounds[i]);
	printf("***************************\n");
	printf("\n");
}



#include "dramrun.c"

int main(int argc, char *argv[])
{
	set_inputs();
	print_inputs();

	char str[12];
	sprintf(str, "%d", input.npar);
	fitfun_initialize(str);

	gsl_rand_init(1);

	double *a = (double *)malloc(input.npar*sizeof(double));
	for (int i = 0; i < input.npar; i++)
		a[i] = 1;

	double Sig[input.npar*input.npar];
	double Lam[input.npar*input.npar];
	covcond(input.c, a, Sig, Lam, input.npar);

	double *mu = (double *)malloc(input.npar*sizeof(double));
	for (int i = 0; i < input.npar; i++)
		mu[i] = 0;

	init_params();
	params.par0    = (double *)malloc(input.npar*sizeof(double));
	for (int i = 0; i < input.npar; i++)
		params.par0[i] = mu[i]+0.1;	// initial value

	params.lbounds = (double *)malloc(input.npar*sizeof(double));
	for (int i = 0; i < input.npar; i++)
		params.lbounds[i] = input.lbounds[i]; //-1e10;	// -inf

	params.ubounds = (double *)malloc(input.npar*sizeof(double));
	for (int i = 0; i < input.npar; i++)
		params.ubounds[i] = input.ubounds[i]; //+1e10;	// -inf

	if (input.pos)
	{
		for (int i = 0; i < input.npar; i++)
			params.lbounds[i] = 0;
	}


	d_mu = mu;
	d_Lam = Lam;

	init_options();
	options.nsimu    = input.nsimu;
	options.adaptint = 500;

	options.qcov = (double *)malloc(input.npar*input.npar*sizeof(double));
	for (int i = 0; i < input.npar; i++)
	for (int j = 0; j < input.npar; j++)
		options.qcov[i*input.npar+j] = Sig[i*input.npar+j]*pow(2.4,2)/input.npar;

	options.drscale  = input.drscale;
	options.adascale = input.adascale; // default is 2.4/sqrt(npar) ;
	options.printint = 2000;

	for (int i = 0; i < input.npar; i++)
		for (int j = 0; j < input.npar; j++)
			printf("options.qcov[%d][%d]=%f\n", i, j, options.qcov[i*input.npar+j]);




	dramrun(input.npar);




	fitfun_finalize();

	return 0;
}




