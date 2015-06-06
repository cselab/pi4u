#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define PROBDIM	2

void fitfun(double const *x, int N, double *res) {
	int i;
	double f = 0.0;   

	for (i=0; i<N-1; i++)   /* rosenbrock */
		f = f + 100.0*pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);

	usleep(100);

	*res = f;
	return; /*return sum;*/
}

int main(int argc, char *argv[])
{
	int rank, size;
	int i;
	double lval = 0, gval;
	double TP[PROBDIM];

	/* print argv numbers */
	for (i = 0; i < argc; i++) {
		printf("simcode: arg %d argv %s\n", i, argv[i]); fflush(0);
	}

	/* read input parameters (argv) into TP */
	//for (i = 1; i < argc; i++) {
	//	TP[i-1] = atof(argv[i]);
	//}

	FILE *fp = fopen("params.dat", "r");
	fscanf(fp, "%lf", &TP[0]);
	fscanf(fp, "%lf", &TP[1]);
	fclose(fp);

	/* run the "simulation" */
	fitfun(TP, PROBDIM, &gval);

	/* write the results in the "fitness" file */
	char fname[256];
	FILE *fd;
	strcpy(fname,"fitness");
	fd = fopen(fname, "w");
	if (fd == NULL) printf("error with fopen\n"); fflush(0);

	//fprintf(fd, "RESULT %.16lf\n", gval);
	fprintf(fd, "%.16lf\n", gval);
	fclose(fd);	

	return 0;
}
