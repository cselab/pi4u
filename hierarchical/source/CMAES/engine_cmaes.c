#include <stdio.h>
#include <stdlib.h> 
#include "cmaes_interface.h"
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "fitfun.h"
#include "../priors/priors.h"

#if defined(_USE_TORC_)

	#include <mpi.h>
	#include <torc.h>

#else

	#include <sys/time.h>
	static double torc_gettime()
	{
        struct timeval t;
        gettimeofday(&t, NULL);
        return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
	}

#endif


#define VERBOSE 1
#define _STEALING_
#define _IODUMP_ 1
#define JOBMAXTIME    0
#define _RESTART_





void set_bounds( double **p_lower_bound, double **p_upper_bound, int dim);
void taskfun(double *x, int *pn, double *res, int *info);
int is_feasible(double *pop, double *lower_bound, double *upper_bound, int dim);

double load_pop_from_file( int step, double * const* pop, double *arFunvals, int dim, int lambda, int *checkp);

void make_all_points_feasible( cmaes_t *evo, double * const *pop, double * lower_bound, double * upper_bound );

double evaluate_population( cmaes_t *evo, double *arFunvals, double * const* pop, Density *d, int step );
void print_the_best( cmaes_t evo, int step );
void write_pop_to_file( cmaes_t evo, double *arFunvals, double * const* pop, int step );

int is_there_enough_time( double gt0, double dt );




int main(int argn, char **args)
{
    cmaes_t evo; 
    double *arFunvals, *const*pop;
    int lambda, dim;
    double gt0, gt1, gt2, gt3;
    double stt = 0.0, dt;
    char dim_str[12];
    int step = 0;
	
	static int checkpoint_restart = 0;

	double *lower_bound, *upper_bound;




	#if defined(_USE_TORC_)
    	torc_register_task(taskfun);
    	torc_init(argn, args, MODE_MS);
	#endif


    gt0 = torc_gettime();


    if(  argn==2  &&  !strcmp(args[1], "-cr")  )	checkpoint_restart = 1;

    
	arFunvals = cmaes_init(&evo, 0, NULL, NULL, 0, 0, "cmaes_initials.par");
    printf("%s\n", cmaes_SayHello(&evo));
    cmaes_ReadSignals(&evo, "cmaes_signals.par");  


    dim = cmaes_Get(&evo, "dim");
  	lambda = cmaes_Get(&evo, "lambda");
	set_bounds( &lower_bound, &upper_bound, dim );


	// Initialize prior distributions
	Density *priors;
	int Nprior;
	read_priors( &priors, &Nprior );

	if( Nprior != dim ){
		printf("The dimension of the prior is different from the dimension of the problem. Exit...");
		exit(1);
	}


	// Initialize log-likelihood
    sprintf(dim_str, "%d", dim );
    fitfun_initialize( dim_str );


    gt1 = torc_gettime();
    while( !cmaes_TestForTermination(&evo) ){

        pop = cmaes_SamplePopulation(&evo); 

        if( checkpoint_restart ){
			dt = load_pop_from_file( step, pop, arFunvals, dim, lambda, &checkpoint_restart);
			stt += dt;
    	}

        
		if( !checkpoint_restart ){
			make_all_points_feasible( &evo, pop, lower_bound, upper_bound );
			dt = evaluate_population( &evo, arFunvals, pop, priors, step );
			stt += dt;
        }

        cmaes_UpdateDistribution(&evo, arFunvals);

        cmaes_ReadSignals(&evo, "cmaes_signals.par"); fflush(stdout);

		print_the_best( evo, step );

		
       	if (!checkpoint_restart){
			write_pop_to_file( evo, arFunvals, pop, step );
        }

		#if defined(_RESTART_)
        	cmaes_WriteToFile(&evo, "resume", "allresumes.dat");
		#endif

		if( ! is_there_enough_time( gt0, dt ) ){
			evo.sp.stopMaxIter=step+1;
			break;
		}
        
		step++;
    }
	gt2 = torc_gettime();

	printf("Stop:\n %s \n",  cmaes_TestForTermination(&evo)); /* print termination reason */
    cmaes_WriteToFile( &evo, "all", "allcmaes.dat" );         /* write final results */
	cmaes_exit(&evo); /* release memory */

    gt3 = torc_gettime();
    
	printf("Total elapsed time = %.3lf  seconds\n", gt3-gt0);
    printf("Initialization time = %.3lf  seconds\n", gt1-gt0);
    printf("Processing time = %.3lf  seconds\n", gt2-gt1);
    printf("Funtion Evaluation time = %.3lf  seconds\n", stt);
    printf("Finalization time = %.3lf  seconds\n", gt3-gt2);




	#if defined(_USE_TORC_)
    	torc_finalize();
	#endif
    
		
	return 0;
}












/*
Assumptions: the feasible domain is convex, the optimum is
     		 not on (or very close to) the domain boundary, initialX is
           	 feasible and initialStandardDeviations are sufficiently small
           	 to prevent quasi-infinite looping. 
*/





//===========================================================================================
//
//
// the function to be minimized
void taskfun(double *x, int *n, double *res, int *info)
{
    (*res) = - fitfun(x, *n, (void *)NULL, info);    // minus for minimization

}






//===========================================================================================
//
//
int is_feasible(double *pop, double *lower_bound, double *upper_bound, int dim)
{
    int i, good;
    for (i = 0; i < dim; i++) {
        good = (lower_bound[i] <= pop[i]) && (pop[i] <= upper_bound[i]);
        if (!good) {
            return 0;
        }
    }
    return 1;
}






//===========================================================================================
//
//
void set_bounds( double **p_lower_bound, double **p_upper_bound, int dim){

    double *lower_bound = malloc(dim*sizeof(double));
    double *upper_bound = malloc(dim*sizeof(double));

    FILE *f = fopen("cmaes_bounds.par", "r");
    
	if (f != NULL){
      
		printf("Reading the bounds from cmaes_bounds.par\n");

      	char line[256];
      	int found;
      	int line_no = 0;
      	for (int i = 0; i < dim; i++) {
        	
			found = 0;
        	while (fgets(line, 256, f)!= NULL) {
          	line_no++;

          	if ((line[0] == '#')||(strlen(line)==0)) continue;

          		char bound[32];
          		sprintf(bound, "B%d", i);
          		if (strstr(line, bound) != NULL) {
            		sscanf(line, "%*s %lf %lf", &lower_bound[i], &upper_bound[i]);
            		found = 1;
            		break;
         		}
        	}
        	if (!found) {
          		printf("Bounds for parameters %d not found in 'cmaes_bounds.par. Exit...'\n", i);
				exit(1);
        	}
        	rewind(f);
        	line_no = 0;
		}
      	fclose(f);
	}
    else {
      printf("Parameters file 'cmaes_bounds.par' could not be opened. Exit...\n");
	  exit(1);
    }


	if( VERBOSE ){
    	printf("Parameter Bounds:\n");
    	for (int i = 0; i < dim; i++) {
			printf("B%d: %15.6f %15.6f\n", i, lower_bound[i], upper_bound[i]);
    	}
	}

	(*p_lower_bound) = lower_bound;
	(*p_upper_bound) = upper_bound;

}







//===========================================================================================
//
//
double load_pop_from_file( int step, double * const* pop, double *arFunvals, int dim, int lambda, int * checkp){
	
 	char filename[256];
	sprintf(filename, "curgen_db_%03d.txt", step);
	FILE *fp = fopen(filename, "r");
    double tt0, tt1 ;
     
  	tt0 = torc_gettime();

	if( fp != NULL ){
   
		for( int i = 0; i < lambda; i++ ){
        	for( int j = 0; j < dim; j++ ){
 				int r = fscanf(fp, "%le", &pop[i][j]);
   
				if(VERBOSE) printf("[%d] pop[%d][%d]=%f\n", r, i, j, pop[i][j]);
                 
				if (r < 0){
					printf("Error occured while reading the (%d,%d) element from %s. Exit...\n",i,j,filename);
					exit(1);
				}
			}
         	
			int r = fscanf(fp, "%le", &arFunvals[i]);
         	
			if(VERBOSE) printf("[%d] arFunvals[%d] = %f\n", r, i, arFunvals[i]);
			
			if (r < 0){
				printf("Error occured while reading the (%d) function value from %s. Exit...\n",i, filename);
				exit(1);
			}
  		}
    	fclose(fp);
 	}
	else
		*checkp = 0;

  	tt1 = torc_gettime();

	return tt1-tt0;

}






//===========================================================================================
//
//
void make_all_points_feasible( cmaes_t *evo, double* const *pop, double * lower_bound, double * upper_bound ){

  	int lambda = cmaes_Get( evo, "lambda");
    int dim    = cmaes_Get( evo, "dim");

	for( int i=0; i<lambda; ++i)
    	while( !is_feasible( pop[i],lower_bound,upper_bound,dim ) )
			cmaes_ReSampleSingle( evo, i );

}






//===========================================================================================
//
//
double evaluate_population( cmaes_t *evo, double *arFunvals, double * const* pop, Density *d, int step ){

  	int lambda = cmaes_Get( evo, "lambda");
    int dim    = cmaes_Get( evo, "dim");
	int info[4];
    double tt0, tt1 ;
    	
  	tt0 = torc_gettime();
	
	for( int i = 0; i < lambda; ++i){
		info[0] = 0; info[1] = 0; info[2] = step; info[3] = i;     /* gen, chain, step, task */
		
		#if defined(_USE_TORC_)
			torc_create( -1, taskfun, 4,
						 dim, MPI_DOUBLE, CALL_BY_VAL,
                       	 1, MPI_INT, CALL_BY_COP,
                		 1, MPI_DOUBLE, CALL_BY_RES,
                    	 4, MPI_INT, CALL_BY_COP,
                 		 pop[i], &dim, &arFunvals[i], info);
		#else
			taskfun(pop[i], &dim, &arFunvals[i], info);
		#endif
	}
	
	#if defined(_USE_TORC_)
		#if defined(_STEALING_)
  			torc_enable_stealing();
		#endif
			torc_waitall();
		#if defined(_STEALING_)
			torc_disable_stealing();
		#endif
	#endif
  	

	// subtruct the log-prior from the log-likelohood
	for( int i=0; i<lambda; i++){
		arFunvals[i] -= prior_log_pdf(d, dim, pop[i]);
	}


	tt1 = torc_gettime();
  
	return tt1-tt0;
}








//===========================================================================================
//
//
void print_the_best( cmaes_t evo, int step ){
	
	#if VERBOSE
		int dim    = cmaes_Get( &evo, "dim");
    	
		const double *xbever = cmaes_GetPtr(&evo, "xbestever");
      	double        fbever = cmaes_Get(   &evo, "fbestever");

      	printf("BEST @ %5d: ", step);
		for( int i = 0; i < dim; i++ )
			printf("%25.16lf ", xbever[i]);
    	printf("%25.16lf\n", fbever);
	#endif
}






//===========================================================================================
//
//
void write_pop_to_file( cmaes_t evo, double *arFunvals, double * const* pop, int step ){

	#if _IODUMP_
		int dim    = cmaes_Get( &evo, "dim");
		int lambda = cmaes_Get( &evo, "lambda");
    	
		char filename[256];
		sprintf(filename, "curgen_db_%03d.txt", step);
		FILE *fp = fopen(filename, "w");
     	for (int i = 0; i < lambda; i++) {
 			for ( int j = 0; j < dim; j++) 
				fprintf(fp, "%.6le ", pop[i][j]);
			fprintf(fp, "%.6le\n", arFunvals[i]);
		}
		fclose(fp);
	#endif

}





//===========================================================================================
//
//
int is_there_enough_time( double gt0, double dt ){
	
	#if (JOBMAXTIME > 0)
	
	double lastgen_time = dt;
	long runt, remt;
    	
	static long maxt = JOBMAXTIME;    //get_job_maxTime(); // from lsf or provided by the user
    printf("job maxtime = %ld\n", maxt);

	runt = torc_gettime()-gt0;    //runt = get_job_runTime(); // from lsf or provided by the application: runt = omp_get_wtime()-gt0;
	remt = maxt - runt;
	printf("job runtime = %ld remtime = %ld\n", runt, remt);

 	if ((lastgen_time*1.1) > remt){
		printf("No more available time, exiting...\n");
		return 0;
	}

	#endif

	return 1;

}





