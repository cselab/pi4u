#include <math.h>
#include <stdio.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>







#if defined(_USE_TORC_)
	#include <mpi.h>

	#ifdef __cplusplus
		extern "C"
		{
	#endif

	#include <torc.h>

	#ifdef __cplusplus
		}
	#endif

#else

	#include <pthread.h>
	static int torc_node_id() { return 0; }
	static int torc_num_nodes() { return 1; }

	#if defined(_USE_OPENMP_)
		#include <omp.h>
		static int torc_i_worker_id() { return omp_get_thread_num(); }
		static int torc_i_num_workers() { return omp_get_max_threads(); }
	#else
		static int torc_i_worker_id() { return 0; }
		static int torc_i_num_workers() { return 1; }
	#endif

#endif



//======================================================================================
//
//


const gsl_rng_type   *T;
gsl_rng   **r;
int   *local_seed;


void gsl_rand_init(int seed){

    int i, local_workers = torc_i_num_workers();
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = (gsl_rng **)malloc(local_workers*sizeof(gsl_rng *));
    for (i = 0; i < local_workers; i++) {
        r[i] = gsl_rng_alloc (T);
    }

    if (seed == 0) seed = time(0);

    for (i = 0; i < local_workers; i++) {
		#if VERBOSE
        printf("node %d: initializing rng %d with seed %d\n", torc_node_id(), i, seed+i+local_workers*torc_node_id());
		#endif
        gsl_rng_set(r[i], seed+i+local_workers*torc_node_id());
    }

    local_seed = (int *)malloc(local_workers*sizeof(int));
    for (i = 0; i < local_workers; i++) {
        local_seed[i] = seed+i+local_workers*torc_node_id();
    }
}



//======================================================================================
//
//
double normal_pdf(double x, double *p){

	return gsl_ran_gaussian_pdf( x-p[0] , p[1] );
}


double normal_log_pdf(double x, double *p){

	double res  =  -0.5*log(2*M_PI)
		   		   -log(p[1]) 
	       		   -0.5 * pow( (x-p[0])/p[1], 2);

	return res;
}


double normal_rnd( double *p ){
    
	double res;
    
	int me = torc_i_worker_id();
    res = p[0] + gsl_ran_gaussian( r[me], p[1] );

    return res;
}


//======================================================================================
//
//
double uniform_pdf(double x, double *p){

	return gsl_ran_flat_pdf( x, p[0] , p[1] );
}


double uniform_log_pdf(double x, double *p){

	if( x>=p[0] && x<=p[1] )
		return -log(p[1]-p[0]);
	else
	   return -INFINITY;	
}


double uniform_rnd( double *p ){
    
	double res;
    
	int me = torc_i_worker_id();
    res = gsl_ran_flat( r[me], p[0], p[1] );

    return res;
}


//======================================================================================
//
//
double exp_pdf(double x, double *p){

	return gsl_ran_exponential_pdf( x, p[0] );
}


double exp_log_pdf(double x, double *p){

	if( x<0 )
		return -INFINITY;	
	else
		return -log(p[0]) - x/p[0];
}


double exp_rnd( double *p ){
    
	double res;
    
	int me = torc_i_worker_id();
    res = gsl_ran_exponential( r[me], p[0] );

    return res;
}



//======================================================================================
//
//      p(x) dx = {1 \over \Gamma(a) b^a} x^{a-1} e^{-x/b} dx
//
double gamma_pdf(double x, double *p){

	return gsl_ran_gamma_pdf( x, p[0], p[1] );
}


double gamma_log_pdf(double x, double *p){

	if( x<0 )
	   return -INFINITY;	

	else{
		double res =  - gsl_sf_lngamma(p[0])
					  - p[0]*log(p[1])
					  + (p[0]-1)*log(x)
					  - x/p[1];

		return res;
	}
}


double gamma_rnd( double *p ){
    
	double res;
    
	int me = torc_i_worker_id();
    res = gsl_ran_gamma( r[me], p[0], p[1] );

    return res;
}










//======================================================================================
//
//    
//

// normal distribution random number N(mu,sigma^2) 
double normalrand(double mu, double sigma){

	int 	me 	= torc_i_worker_id();
    double	res = gsl_ran_gaussian( r[me], sigma );

	return res;
}



// uniform distribution random number between a and b
double uniformrand(double a, double b){

	int 	me 	= torc_i_worker_id();
    double	res = gsl_ran_flat( r[me], a, b );

	return res;
}



