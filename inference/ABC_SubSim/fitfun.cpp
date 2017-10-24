#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#if defined(_USE_TORC_)
#include <mpi.h>
extern "C"
{
#include <torc.h>
}
#else
int torc_worker_id()
{
	return 0;
}
#endif

static int dbg_display = 0;

double fitfun(double *x, int n, int info[4])
{
	int i;

	char taskname[256];
	//sprintf(taskname, "tmpdir.%d.%d", getpid(), torc_worker_id());
	sprintf(taskname, "tmpdir.%d.%d.%d.%d", info[0], info[1], info[2], info[3]);

	if (dbg_display) {
		printf("worker(%d): running task %s with params (", torc_worker_id(), taskname);
		for (i = 0; i < n-1; i++)
			printf("%.4lf,", x[i]);
		printf("%.4lf)\n", x[i]);
		fflush(0);
	}

	double f = 0.0;
	for (i=0; i<n-1; i++)   /* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
//	f = -f;

	double res = log(f);
//	double res = f;
	return res;
}
