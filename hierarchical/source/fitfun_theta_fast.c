#include "fitfun.h"


#define DATABASE  "."
#define DATAFILE "data.txt"

#define BUFLEN 1024

#define NREC 10000





// number of data in data file
static int Nd=0;
static double *datax, *datay;


void fitfun_initialize() {

	FILE *fp;

	// open the data file
	char filename[BUFLEN];
	snprintf(filename, BUFLEN, "%s/%s", DATABASE, DATAFILE);

  	fp = fopen(filename, "r");
	if (fp == NULL) {
   		printf("\n %s  does not exist. Exiting...\n", filename);
  		exit(1);
   	}
	printf("\nData file %s succesfully opened \n",filename);

	// read the number of lines of data file
	char ch;

 	while (!feof(fp)) {
   		ch = fgetc(fp);
   		if (ch == '\n')  Nd++;
  	}
  	rewind(fp);


	datax = (double*)malloc(Nd*sizeof(double));
	datay = (double*)malloc(Nd*sizeof(double));
	for(int i=0; i<Nd; i++){
		fscanf( fp, "%lf", &datax[i]);
		fscanf( fp, "%lf", &datay[i]);
	}

	printf("\n%d data succesfully read from %s \n",Nd,filename);
	
	
	fclose(fp);

}



void fitfun_finalize() {
    
	free(datax);
	free(datay);

}


void my_model(double *x, double *y, double *c, double n){

	for(int i=0; i<n; i++)
		y[i] = c[0] * sin( c[1]*x[i]+c[2] );

}




double fitfun(double *x, int n, void *output, int *info) {


	
	/*
	x[0]=2.;
	x[1]=3.;
	x[2]=1.;
	x[3]=0.5;
	*/

    double res;
    double sigma2 = pow(x[n-1],2);

	double *y;
	y = (double*)malloc(Nd*sizeof(double));
	my_model(datax,y,x,Nd);


	double ssn=0.;
	for(int i=0; i<Nd; i++)
		ssn += pow( datay[i] - y[i], 2);

	res = Nd*log(2*M_PI) + Nd*log(sigma2) + ssn/sigma2;
	res *= -0.5;

	/*
	printf("\n Nd = %d \n",Nd);
	printf("\n datay(1) = %lf \n",datay[0]);
	printf("\n datay(Nd) = %lf \n",datay[Nd-1]);
	printf("\n y(1) = %lf \n",y[0]);
	printf("\n y(Nd) = %lf \n",y[Nd-1]);
	printf("\n %lf \n",res);
	exit(1);
	*/

    return res;
}



