
#ifndef _TMCMC_STATS_H_
#define _TMCMC_STATS_H_

void calculate_statistics(double flc[], unsigned int n, int nselections, int gen, unsigned int sel[]);


typedef struct fparam_s {
	double *fj;
	int     fn;
	double  pj;
	double  tol;
} fparam_t;



#endif
