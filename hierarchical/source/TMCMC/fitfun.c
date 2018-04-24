#include <math.h>


void fitfun_initialize(char *s){
}



void fitfun_finalize(){
}




double fitfun(double *x, int N, void *output, int *info){
	double f;

	int i;
	f = 0.0;
	for (i=0; i<N-1; i++)	/* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
	f = -f;

	return f;
}
