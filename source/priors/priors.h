#ifndef _PRIORS_H_
#define _PRIORS_H_



typedef struct{
	double (*f)(double, double *);		// Density function
	double (*lf)(double, double *);		// Log of density
	double (*r)( double *);				// Random number
	double *par;						// Parameters of density
	int npar;							// Number of parameters
	char   name[24];					// Name of the distribution
} Density;


//TODO: define the struct Prior:
//typedef struct{
//	Density *d;
//	int Nd;
//} Prior;


void delete_density( Density *d );
void delete_prior( Density *d, int N);

double eval_density(Density d, double x);
double eval_log_density(Density d, double x);
double eval_random( Density d );

void  print_density( Density d );
void  print_priors( Density *d, int N);

double prior_pdf( Density *d, int N, double *x);
double prior_log_pdf( Density *d, int N, double *x);

void read_priors(const char *file,  Density **p_priors, int *p_N );

void reassign_prior( Density *p, int Np, double *psi );
void new_prior_from_prior( Density **new_prior, Density *from_prior, int Npr );

void check_n( int N );



#endif
