#include <math.h>
#include <unistd.h>

void fitfun(double *x, int nx, void *output, int *winfo, double *result, int ncon, double *constraints) // ny must be also available
{
	usleep(10*1000);

//	result[0] = pow((x[0]+2.0),2.0)-10;
//	result[1] = pow((x[0]-2.0),2.0)+20;

	double y = 0;

	int i;
	for (i = 0; i < nx; i++) y += x[i];

	result[0] = pow((y+2.0),2.0)-10;
	result[1] = pow((y-2.0),2.0)+20;
}

