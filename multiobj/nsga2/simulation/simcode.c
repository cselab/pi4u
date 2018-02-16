#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#define PROBDIM	1
#define NOBJ	2
#define NCONSTR	1

#define PARAMSFILE      "params.txt"
#define OUTPUTFILE      "fitness.txt"
#define CONSTRFILE      "constr_viol.txt"

void fitfun(double *x, int *pnx, double *obj, int *pno, double *constr, int *pncon)
{
	int nx = *pnx;
	int no = *pno;
	int ncon = *pncon;

	double y = 0;
	for (int i = 0; i < nx; i++) y += x[i];

	obj[0] = pow((y+2.0),2.0)-10;
	obj[1] = pow((y-2.0),2.0)+20;

	constr[0] = obj[1] - 25;

	return;
}

int main(int argc, char *argv[])
{
	int nx = PROBDIM;
	double x[PROBDIM];
	int no = NOBJ;
	double obj[NOBJ];
	int ncon = NCONSTR;
	double constr[NCONSTR+1];

	/* print argv numbers */
	for (int i = 0; i < argc; i++) {
		printf("simcode: arg %d argv %s\n", i, argv[i]); fflush(0);
	}

	/* read input parameters (argv) into TP */

	FILE *fp = fopen(PARAMSFILE, "r");
	if (fp == NULL) {
		printf("Missing %s file. Exiting...\n", PARAMSFILE);
		exit(1);
	}

	for (int i = 0; i < nx; i++)
		fscanf(fp, "%lf", &x[i]);
	fclose(fp);

	/* run the "simulation" */
	for (int i = 0; i < nx; i++) {
		printf("input parameter %d = %f\n", i, x[i]);
	}

	fitfun (x, &nx, obj, &no, constr, &ncon);

	/* write the results in the "OUTPUT" file */
	char foname[256];
	FILE *fo;
	strcpy(foname,OUTPUTFILE);
	fo = fopen(foname, "w");
	if (fo == NULL) {
		printf("Could not create  %s file. Exiting...\n", OUTPUTFILE);
		exit(1);
	}

	for (int i = 0; i < no; i++)
		fprintf(fo, "%.16lf\n", obj[i]);
	fclose(fo);


	/* write the constraints violations in the "CONSTRFILE" file */
	char fcname[256];
	FILE *fc;
	strcpy(fcname,CONSTRFILE);
	fc = fopen(fcname, "w");
	if (fc == NULL) {
		printf("Could not create  %s file. Exiting...\n", CONSTRFILE);
		exit(1);
	}

	for (int i = 0; i < ncon; i++)
		fprintf(fc, "%.16lf\n", constr[i]);
	fclose(fc);	


	return 0;
}
