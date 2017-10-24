#include <math.h>

// activate one of the following options
//#define _USE_ROSENBROCK_
//#define _USE_BVNPDF_
#define _USE_MIXED_BVNPDF_
//#define _USE_MIXED_MVNPDF_

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
	return P;
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

double fitfun(double /*const*/ *x, int N, void *output, int *info)
{
	double f;

#if defined(_USE_ROSENBROCK_)
	int i;
	f = 0.0;
	for (i=0; i<N-1; i++)	/* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
	f = -f;
//	f = -log(f);	// peh xxx: logval = 1
#endif

#if defined(_USE_BVNPDF_)
	f = log(bvnpdf(x, N));
#endif

#if defined(_USE_MIXED_BVNPDF_)
	f = log(mixedbvnpdf(x, N));
#endif

#if defined(_USE_MIXED_MVNPDF_)
	f = log(mixedmvnpdf(x, N));
#endif
	return f;
}

