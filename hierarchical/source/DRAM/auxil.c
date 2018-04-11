#include "auxil.h"

#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>





const gsl_rng_type *T;
gsl_rng *r;
static pthread_mutex_t _rm = PTHREAD_MUTEX_INITIALIZER;

void gsl_rand_init(int seed)
{
		gsl_rng_env_setup();

		T = gsl_rng_default;
		r = gsl_rng_alloc (T);
		if (seed == 0) 
			gsl_rng_set(r, time(0) );
		else
			gsl_rng_set(r, seed);
}



/* normal distribution random number N(mu,rho^2)*/
double normalrand(double mu, double var)
{
	double res;

	pthread_mutex_lock(&_rm);
	res = mu + gsl_ran_gaussian(r, var);
	pthread_mutex_unlock(&_rm);
	return res;

/*	return mu + gsl_ran_gaussian(r, var);*/
}



/* uniform (flat) distribution random number between a and b */
double uniformrand(double a, double b)
{
	double res;

	pthread_mutex_lock(&_rm);
	res = gsl_ran_flat(r, a, b);
	pthread_mutex_unlock(&_rm);
	return res;
}



