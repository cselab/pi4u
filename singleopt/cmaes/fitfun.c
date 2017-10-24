#include <math.h>

// activate one of the following options
#define _USE_ROSENBROCK_
//#define _USE_RASTRIGIN_
//#define _USE_CIGTAB_

double fitfun(double /*const*/ *x, int N, void *output, int *info)
{
	int i;
	double f;

#if TESTING
	printf("fitfun {%d,%d,%d,%d}\n", winfo[0], winfo[1], winfo[2], winfo[3]);
	usleep(100*1000);
#endif


#if defined(_USE_ROSENBROCK_)
	f = 0.0;
	for (i=0; i<N-1; i++)	/* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
	f = -f;
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

	return f;
}

