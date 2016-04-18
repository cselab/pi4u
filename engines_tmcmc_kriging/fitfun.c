#include <math.h>

// activate one of the following options
//#define _USE_CIRCLE_
//#define _USE_ROSENBROCK_
//#define _USE_RASTRIGIN_
//#define _USE_CIGTAB_
//#define _USE_BVNPDF_
//#define _USE_SUBSETFUN_
//#define _USE_HIMMELBLAU_
//#define _USE_PA_DEMO_
#define _USE_LOGNORMPDF_
//#define _USE_MIXED_BVNPDF_
//#define _USE_MIXED_MVNPDF_
//#define _USE_MVNPDF_

#if defined(_USE_MIXED_BVNPDF_)
#include "gsl_headers.h"
double mixedbvnpdf(double *x, int n) /* bivariate */
{
	double P;

	P = gsl_ran_bivariate_gaussian_pdf(x[0]-(-5), x[1]-(-5), 1, 1, 0);
	P += gsl_ran_bivariate_gaussian_pdf(x[0]-(+5), x[1]-(+5), 1, 1, 0);
	return P;
}
#endif

#if defined(_USE_MIXED_MVNPDF_)
#include "gsl_headers.h"
double mvnpdf(int n, double *xv, double *mv, double *vm);
double mixedmvnpdf(double *x, int n) /* multivariate */
{
	double P = 0;
	double m1[n];
	double m2[n];

	int i;
	for (i = 0; i < n; i++) m1[i] = -5;
	for (i = 0; i < n; i++) m2[i] = +5;

	P =  mvnpdf(n, x, m1, NULL);
	P += mvnpdf(n, x, m2, NULL);
	return log(P);
}
#endif

#if defined(_USE_BVNPDF_)
#include "gsl_headers.h"
double bvnpdf(double *x, int n) /* bivariate */
{
	double P;

	P = gsl_ran_bivariate_gaussian_pdf(x[0], x[1], 1, 1, 0);
	return P;
}

#endif

#if defined(_USE_MVNPDF_)
#include "gsl_headers.h"
double mvnpdf(int n, double *xv, double *mv, double *vm);
#endif


double fitfun(double /*const*/ *x, int N, void *output, int *info)
{
	double f;

//	printf("fitfun {%d,%d,%d,%d}\n", winfo[0], winfo[1], winfo[2], winfo[3]);
//	usleep(100*1000);


#if defined(_USE_CIRCLE_)
	int i;
	f = 0.0;
	for(i = 0; i < N; i++) 
		f += x[i]*x[i]; 
	f = -f;
#endif


#if defined(_USE_ROSENBROCK_)
	int i;
	f = 0.0;
	for (i=0; i<N-1; i++)	/* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
	f = -f;
//	f = -log(f);	// peh xxx: logval = 1
#endif

#if defined(_USE_RASTRIGIN_)
	int i;
	const double pi=3.14159265358979;
	f = 0.0;
	for (i=0; i<N; i++)	/* rastrigin */
		f = f + pow(x[i],2) + 10.0 - 10.0*cos(2*pi*x[i]);
	f = -f;
#endif

#if defined(_USE_CIGTAB_)
	f = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1];	/* cigtab */
	for(i = 2; i < N; i++) 
		f += x[i]*x[i]; 
	f= -f;
#endif


#if defined(_USE_BVNPDF_)
	f = log(bvnpdf(x, N));
#endif

#if defined(_USE_SUBSETFUN_)
	f = 8*exp(-(pow(x[0],2.0)+pow(x[1],2.0))) + 2*exp(-(pow(x[0]-5,2.0)+pow(x[1]-4,2.0))) + 1.0 + (x[0]*x[1])/10.0;
#endif

#if defined(_USE_HIMMELBLAU_)
//	f = (x^2+y-11)^2 + (x+y^2-7)^2
//	f = 	pow(x[0]*x[0]+x[1]-11, 2) + 
//		pow(x[0]+x[1]*x[1]- 7, 2);
//	f = -0.1*f;

	int i;
	//evaluate Himmel-blau function at given point x
	f = 0.0;
	for(i=0; i<N-1; i++) {
		//f(x,y) =       (x^2 + y - 11)^2      +   (x + y^2 - 7.0)^2
		// /  |  |
		f = f + pow((x[i]*x[i]+x[i+1]-11.0),2) + pow((x[i]+x[i+1]*x[i+1]-7.0),2);
	}
	f = (-0.1) * f; //exp((-0.1)*f); //do not use exp as you assumed that finally you used log(exp(f)) = f
#endif

#if defined(_USE_PA_DEMO_)
//	PA & GT
	int i;
	if (N != 2) { printf("N must be 2! Aborting execution\n"); exit(1); } // assert(N==2);
	double x0 = x[0];
	double alpha = x[1];
	double t[] = {0.0600, 10.0100, 20.0200, 30.0300, 40.0400, 50.0501};
	double data[] = { 1.6507, -0.0647, -0.8515, -0.1902, -0.6457, 0.0340};
	double sigma2 = 0.5;
	int Ndata = sizeof(t)/sizeof(t[0]);
	f = 0;
	for (i = 0; i < Ndata; ++i) {
		const double ft = data[i] - x0*exp(-alpha*t[i]);
		f -= ft*ft;
	}
	f = 0.5*f/sigma2;
#endif

#if defined(_USE_NORMPDF_)
	int i;
	//evaluate the nxn normal distribution
	const double pi=3.14159265358979;
	f = 0.0;
	for(i = 0; i < N; i++)
		f = f + pow(x[i],2.0);
	
	f = 1.0/pow(2*pi,0.5*N) * exp(-0.5*f);
#endif


#if defined(_USE_LOGNORMPDF_)
	int i;
	//evaluate the log of a nxn normal distribution
	const double pi=3.14159265358979;
	//double half_N = N / 2.0;
	f = 0.0;
	for(i = 0; i < N; i++){
		f = f + pow(x[i],2.0);
	}
	
	//  (-1/2) * (1/(2pi)^(N/2) * sqrt(det(Covariance))) * (x-mu)^T * Covariance^(-1) * (x-mu)
	f = (-1.0/2.0) * N * log(2*pi) - 0.5*f;
#endif

#if defined(_USE_MIXED_BVNPDF_)
	f = log(mixedbvnpdf(x, N));
//	f = mixedbvnpdf(x, N);
#endif

#if defined(_USE_MIXED_MVNPDF_)
	f = mixedmvnpdf(x, N);
#endif

#if defined(_USE_MVNPDF_)
	{
		double m[N];
	
		int i;
		for (i = 0; i < N; i++) m[i] = 0;

		f = mvnpdf(N, x, m, NULL);
	}
#endif

	return f;
}

