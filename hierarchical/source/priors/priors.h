#ifndef _PRIORS_H_
#define _PRIORS_H_



typedef struct{
	double (*f)(double, double *);		// Density function
	double (*lf)(double, double *);		// Log of density
	double (*r)( double *);				// Random number
	double *par;						// Parameters of density
	char   name[24];					// Name of the distribution
} Density;



double eval_density(Density d, double x);
double eval_log_density(Density d, double x);
double eval_random( Density d );

double print_density( Density d );
double print_priors( Density *d, int N);

double prior_pdf( Density *d, int N, double x);
double prior_log_pdf( Density *d, int N, double x);

void read_priors( Density **p_priors, int *p_N );

void check_n( int N );



#endif
