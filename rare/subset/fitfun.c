#include <math.h>

//#define _USE_ROSENBROCK_
//#define _USE_RASTRIGIN_
#define _USE_SUBSETFUN_

double fitfun(double /*const*/ *x, int N, void *output, int *info)
{
	double f;

#if defined(_USE_ROSENBROCK_)
	int i;
	f = 0.0;
	for (i=0; i<N-1; i++)   /* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
#endif

#if defined(_USE_RASTRIGIN_)
	int i;
	const double pi=3.14159265358979;
	f = 0.0;
	for (i=0; i<N; i++)     /* rastrigin */
		f = f + pow(x[i],2) + 10.0 - 10.0*cos(2*pi*x[i]);
#endif

#if defined(_USE_SUBSETFUN_)
	f = 8*exp(-(pow(x[0],2.0)+pow(x[1],2.0))) + 2*exp(-(pow(x[0]-5,2.0)+pow(x[1]-4,2.0))) + 1.0 + (x[0]*x[1])/10.0;
#endif

	return f;
}

