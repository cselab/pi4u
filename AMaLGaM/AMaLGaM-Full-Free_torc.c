/**
 * AMaLGaM-Full-Free.c
 *
 * Copyright (c) 1998-2010 Peter A.N. Bosman
 *
 * The software in this file is the proprietary information of
 * Peter A.N. Bosman.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 *
 * Adapted Maximum-Likelihood Gaussian Model
 * Iterated Density-Estimation Evolutionary Algorithm
 * with a full covariance matrix and free of parameters
 *
 * In this implementation, minimization is assumed.
 *
 * The software in this file is the result of (ongoing) scientific research.
 * The following people have been actively involved in this research over
 * the years:
 * - Peter A.N. Bosman
 * - Dirk Thierens
 * - Jörn Grahl
 *
 * This is the most up-to-date literature reference regarding this software:
 *
 * P.A.N. Bosman. On Empirical Memory Design, Faster Selection of Bayesian
 * Factorizations and Parameter-Free Gaussian EDAs. In G. Raidl, E. Alba,
 * J. Bacardit, C. Bates Congdon, H.-G. Beyer, M. Birattari, C. Blum,
 * P.A.N. Bosman, D. Corne, C. Cotta, M. Di Penta, B. Doerr, R. Drechsler,
 * M. Ebner, J. Grahl, T. Jansen, J. Knowles, T. Lenaerts, M. Middendorf,
 * J.F. Miller, M. O'Neill, R. Poli, G. Squillero, K. Stanley, T. Stützle
 * and J. van Hemert, editors, Proceedings of the Genetic and Evolutionary
 * Computation Conference - GECCO-2009, pages 389-396, ACM Press, New York,
 * New York, 2009. 
 */

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <torc.h>
#include "fitfun.c"



/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
void *Malloc( long size );
double **matrixNew( int n, int m );
double vectorDotProduct( double *vector0, double *vector1, int n0 );
double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 );
double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 );
int blasDSWAP( int n, double *dx, int incx, double *dy, int incy );
int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy);
void blasDSCAL( int n, double sa, double x[], int incx );
int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] );
double **choleskyDecomposition( double **matrix, int n );
int linpackDTRDI( double t[], int ldt, int n );
double **matrixLowerTriangularInverse( double **matrix, int n );
void matrixWriteToFile( FILE *file, double **matrix, int n0, int n1 );
int *mergeSort( double *array, int array_size );
void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q );
void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q );
int *mergeSortFitness( double *objectives, double *constraints, int number_of_solutions );
void mergeSortFitnessWithinBounds( double *objectives, double *constraints, int *sorted, int *tosort, int p, int q );
void mergeSortFitnessMerge( double *objectives, double *constraints, int *sorted, int *tosort, int p, int r, int q );
void interpretCommandLine( int argc, char **argv );
void parseCommandLine( int argc, char **argv );
void parseOptions( int argc, char **argv, int *index );
void printAllInstalledProblems( void );
void optionError( char **argv, int index );
void parseParameters( int argc, char **argv, int *index );
void printUsage( void );
void checkOptions( void );
void printVerboseOverview( void );
double randomRealUniform01( void );
int randomInt( int maximum );
double random1DNormalUnit( void );
char *installedProblemName( int index );
int numberOfInstalledProblems( void );
double installedProblemLowerRangeBound( int index, int dimension );
double installedProblemUpperRangeBound( int index, int dimension );
short isParameterInRangeBounds( double parameter, int dimension );
void installedProblemEvaluation( int *pindex, double *parameters, double *objective_value, double *constraint_value, int *info );
void installedProblemEvaluationWithoutRotation( int index, double *parameters, double *objective_value, double *constraint_value, int *info );
void rosenbrockFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int *info );
double rosenbrockFunctionLowerRangeBound( int dimension );
double rosenbrockFunctionUpperRangeBound( int dimension );
void initialize( void );
void initializeMemory( void );
void initializeRandomNumberGenerator( void );
void initializeParameterRangeBounds( void );
void initializeDistributionMultipliers( void );
void initializePopulationsAndFitnessValues( void );
void initializeObjectiveRotationMatrix( void );
void computeRanks( void );
void computeRanksForOnePopulation( int population_index );
double distanceInParameterSpace( double *solution_a, double *solution_b );
void writeGenerationalStatistics( void );
void writeGenerationalSolutions( short final );
void writeGenerationalSolutionsBest( short final );
void writeGenerationalSolutionsBestFinal( void );
void writeGenerationalSolutionsBestNotFinal( void );
short checkTerminationCondition( void );
short checkTerminationConditionForRunOnce( void );
short checkNumberOfEvaluationsTerminationCondition( void );
short checkVTRTerminationCondition( void );
void determineBestSolutionInCurrentPopulations( int *population_of_best, int *index_of_best );
void checkFitnessVarianceTermination( void );
short checkFitnessVarianceTerminationSinglePopulation( int population_index );
void checkDistributionMultiplierTerminationCondition( void );
void makeSelections( void );
void makeSelectionsForOnePopulation( int population_index );
void makeSelectionsForOnePopulationUsingDiversityOnRank0( int population_index );
void makePopulations( void );
void estimateParametersAllPopulations( void );
void estimateParameters( int population_index );
void estimateParametersML( int population_index );
void estimateMeanVectorML( int population_index );
void estimateCovarianceMatrixML( int population_index );
void copyBestSolutionsToPopulations( void );
void applyDistributionMultipliers( void );
void generateAndEvaluateNewSolutionsToFillPopulations( void );
void computeParametersForSampling( int population_index );
double *generateNewSolution( int population_index );
void adaptDistributionMultipliers( void );
short generationalImprovementForOnePopulation( int population_index, double *st_dev_ratio );
short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
double getStDevRatio( int population_index, double *parameters );
void determineBestSolutionSoFar( void );
void ezilaitini( void );
void ezilaitiniMemory( void );
void ezilaitiniDistributionMultipliers( void );
void ezilaitiniObjectiveRotationMatrix( void );
void run( void );
void runOnce( void );
int main( int argc, char **argv );
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
short      write_generational_statistics,    /* Whether to compute and write statistics every generation (0 = no). */
           write_generational_solutions,     /* Whether to write the population every generation (0 = no). */
           print_verbose_overview,           /* Whether to print a overview of settings (0 = no). */
           use_vtr,                          /* Whether to terminate at the value-to-reach (VTR) (0 = no). */
           haveNextNextGaussian = 0,         /* Internally used variable for sampling the normal distribution. */
          *populations_terminated;           /* Which populations have terminated (array). */
int        number_of_parameters,             /* The number of parameters to be optimized. */
           population_size,                  /* The size of each population. */
           selection_size,                   /* The size of the selection for each population. */
           maximum_number_of_evaluations,    /* The maximum number of evaluations. */
           number_of_evaluations,            /* The current number of times a function evaluation was performed. */
           number_of_generations,            /* The current generation count. */
           number_of_populations,            /* The number of parallel populations that initially partition the search space. */
           number_of_starts,                 /* The number of times the algorithm was started. */
           problem_index,                    /* The index of the optimization problem. */
          *samples_drawn_from_normal,        /* The number of samples drawn from the i-th normal in the last generation. */
          *out_of_bounds_draws,              /* The number of draws that resulted in an out-of-bounds sample. */
          *no_improvement_stretch,           /* The number of subsequent generations without an improvement while the distribution multiplier is <= 1.0, for each population separately. */
           maximum_no_improvement_stretch,   /* The maximum number of subsequent generations without an improvement while the distribution multiplier is <= 1.0. */
           maximum_number_of_populations;    /* The maximum number of parallel populations that initially partition the search space. */
double     tau,                              /* The selection truncation percentile (in [1/population_size,1]). */
           alpha_AMS,                        /* The percentile of offspring to apply AMS (anticipated mean shift) to. */
           delta_AMS,                        /* The adaptation length for AMS (anticipated mean shift). */
        ***populations,                      /* The populations containing the solutions. */
         **objective_values,                 /* Objective values for population members. */
         **constraint_values,                /* Sum of all constraint violations for population members. */
         **ranks,                            /* Ranks of population members. */
        ***selections,                       /* Selected solutions, one for each population. */
         **objective_values_selections,      /* Objective values of selected solutions. */
         **constraint_values_selections,     /* Sum of all constraint violations of selected solutions. */
          *lower_range_bounds,               /* The respected lower bounds on parameters. */
          *upper_range_bounds,               /* The respected upper bounds on parameters. */
          *lower_init_ranges,                /* The initialization range lower bound. */
          *upper_init_ranges,                /* The initialization range upper bound */
           lower_user_range,                 /* The initial lower range-bound indicated by the user (same for all dimensions). */
           upper_user_range,                 /* The initial upper range-bound indicated by the user (same for all dimensions). */
           rotation_angle,                   /* The angle of rotation to be applied to the problem. */
          *distribution_multipliers,         /* Distribution multipliers (AVS mechanism), one for each population. */
           distribution_multiplier_increase, /* The multiplicative distribution multiplier increase. */
           distribution_multiplier_decrease, /* The multiplicative distribution multiplier decrease. */
           st_dev_ratio_threshold,           /* The maximum ratio of the distance of the average improvement to the mean compared to the distance of one standard deviation before triggering AVS (SDR mechanism). */
           vtr,                              /* The value-to-reach (function value of best solution that is feasible). */
           fitness_variance_tolerance,       /* The minimum fitness variance level that is allowed. */
         **mean_vectors,                     /* The mean vectors, one for each population. */
         **mean_vectors_previous,            /* The mean vectors of the previous generation, one for each population. */
        ***covariance_matrices,              /* The covariance matrices to be used for sampling, one for each population. */
        ***cholesky_factors_lower_triangle,  /* The unique lower triangular matrix of the Cholesky factorization. */
           nextNextGaussian,                 /* Internally used variable for sampling the normal distribution. */
         **rotation_matrix,                  /* The rotation matrix to be applied before evaluating. */
          *best_so_far_solution,             /* Best solution found so far over all starts. */
           best_so_far_objective_value,      /* Objective value of best solution found so far over all starts. */
           best_so_far_constraint_value;     /* Sum-of-constraints value of best solution found so far over all starts. */
int64_t    random_seed,                      /* The seed used for the random-number generator. */
           random_seed_changing;             /* Internally used variable for randomly setting a random seed. */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/




/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Constants -=-=-=-=-=-=-=-=-=-=-=-=-=-*/
#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/






/*-=-=-=-=-=-=-=-=-=-=-= Section Elementary Operations -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Allocates memory and exits the program in case of a memory allocation failure.
 */
void *Malloc( long size )
{
  void *result;

  result = (void *) malloc( size );
  if( !result )
  {
    printf("\n");
    printf("Error while allocating memory in Malloc( %ld ), aborting program.", size);
    printf("\n");

    exit( 0 );
  }

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Matrix -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Creates a new matrix with dimensions n x m.
 */
double **matrixNew( int n, int m )
{
  int      i;
  double **result;

  result = (double **) malloc( n*( sizeof( double * ) ) );
  for( i = 0; i < n; i++ )
    result[i] = (double *) malloc( m*( sizeof( double ) ) );

  return( result );
}

/**
 * Computes the dot product of two vectors of the same dimensionality n0.
 */
double vectorDotProduct( double *vector0, double *vector1, int n0 )
{
  int    i;
  double result;
  
  result = 0.0;
  for( i = 0; i < n0; i++ )
    result += vector0[i]*vector1[i];
    
  return( result );
}

/**
 * Computes the multiplication Av of a matrix A and a vector v
 * where matrix A has dimensions n0 x n1 and vector v has
 * dimensionality n1.
 */
double *matrixVectorMultiplication( double **matrix, double *vector, int n0, int n1 )
{
  int     i;
  double *result;
  
  result = (double *) malloc( n0*sizeof( double ) );
  for( i = 0; i < n0; i++ )
    result[i] = vectorDotProduct( matrix[i], vector, n1 );
  
  return( result );
}

/**
 * Computes the matrix multiplication of two matrices A and B
 * of dimensions A: n0 x n1 and B: n1 x n2.
 */
double **matrixMatrixMultiplication( double **matrix0, double **matrix1, int n0, int n1, int n2 )
{
  int     i, j, k;
  double **result;
  
  result = (double **) malloc( n0*sizeof( double * ) );
  for( i = 0; i < n0; i++ )
    result[i] = (double *) malloc( n2*sizeof( double ) );

  for( i = 0; i < n0; i++ )
  {
    for( j = 0; j < n2; j++ )
    {
      result[i][j] = 0;
      for( k = 0; k < n1; k++ )
        result[i][j] += matrix0[i][k]*matrix1[k][j];
    }
  }
  
  return( result );
}

/**
 * BLAS subroutine.
 */
int blasDSWAP( int n, double *dx, int incx, double *dy, int incy )
{
  double dtmp;
  
  if (n > 0)
  {
    incx *= sizeof( double );
    incy *= sizeof( double );

    dtmp  = (*dx);
    *dx   = (*dy);
    *dy   = dtmp;

    while( (--n) > 0 )
    {
      dx = (double *) ((char *) dx + incx);
      dy = (double *) ((char *) dy + incy);
      dtmp = (*dx); *dx = (*dy); *dy = dtmp;
    }
  }

  return( 0 );
}

/**
 * BLAS subroutine.
 */
int blasDAXPY(int n, double da, double *dx, int incx, double *dy, int incy)
{
  double dtmp0, dtmp, *dx0, *dy0;
  
  if( n > 0 && da != 0. )
  {
    incx *= sizeof(double);
    incy *= sizeof(double);
    *dy  += da * (*dx);

    if( (n & 1) == 0 )
    {
      dx   = (double *) ((char *) dx + incx);
      dy   = (double *) ((char *) dy + incy);
      *dy += da * (*dx);
      --n;
    }
    n = n >> 1;
    while( n > 0 )
    {
      dy0   = (double *) ((char *) dy + incy);
      dy    = (double *) ((char *) dy0 + incy);
      dtmp0 = (*dy0);
      dtmp  = (*dy); 
      dx0   = (double *) ((char *) dx + incx); 
      dx    = (double *) ((char *) dx0 + incx);
      *dy0  = dtmp0 + da * (*dx0);
      *dy   = dtmp + da * (*dx);
      --n;
    }
  }

  return( 0 );
}

/**
 * BLAS subroutine.
 */
void blasDSCAL( int n, double sa, double x[], int incx )
{
  int i, ix, m;

  if( n <= 0 )
  {
  }
  else if( incx == 1 )
  {
    m = n % 5;

    for( i = 0; i < m; i++ )
    {
      x[i] = sa * x[i];
    }

    for( i = m; i < n; i = i + 5 )
    {
      x[i]   = sa * x[i];
      x[i+1] = sa * x[i+1];
      x[i+2] = sa * x[i+2];
      x[i+3] = sa * x[i+3];
      x[i+4] = sa * x[i+4];
    }
  }
  else
  {
    if( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    for( i = 0; i < n; i++ )
    {
      x[ix] = sa * x[ix];
      ix = ix + incx;
    }
  }
}

/**
 * LINPACK subroutine.
 */
int linpackDCHDC( double a[], int lda, int p, double work[], int ipvt[] )
{
  int    info, j, jp, k, l, maxl, pl, pu;
  double maxdia, temp;

  pl   = 1;
  pu   = 0;
  info = p;
  for( k = 1; k <= p; k++ )
  {
    maxdia = a[k-1+(k-1)*lda];
    maxl   = k;
    if( pl <= k && k < pu )
    {
      for( l = k+1; l <= pu; l++ )
      {
        if( maxdia < a[l-1+(l-1)*lda] )
        {
          maxdia = a[l-1+(l-1)*lda];
          maxl   = l;
        }
      }
    }

    if( maxdia <= 0.0 )
    {
      info = k - 1;

      return( info );
    }

    if( k != maxl )
    {
      blasDSWAP( k-1, a+0+(k-1)*lda, 1, a+0+(maxl-1)*lda, 1 );

      a[maxl-1+(maxl-1)*lda] = a[k-1+(k-1)*lda];
      a[k-1+(k-1)*lda]       = maxdia;
      jp                     = ipvt[maxl-1];
      ipvt[maxl-1]           = ipvt[k-1];
      ipvt[k-1]              = jp;
    }
    work[k-1]        = sqrt( a[k-1+(k-1)*lda] );
    a[k-1+(k-1)*lda] = work[k-1];

    for( j = k+1; j <= p; j++ )
    {
      if( k != maxl )
      {
        if( j < maxl )
        {
          temp                = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda]    = a[j-1+(maxl-1)*lda];
          a[j-1+(maxl-1)*lda] = temp;
        }
        else if ( maxl < j )
        {
          temp                = a[k-1+(j-1)*lda];
          a[k-1+(j-1)*lda]    = a[maxl-1+(j-1)*lda];
          a[maxl-1+(j-1)*lda] = temp;
        }
      }
      a[k-1+(j-1)*lda] = a[k-1+(j-1)*lda] / work[k-1];
      work[j-1]        = a[k-1+(j-1)*lda];
      temp             = -a[k-1+(j-1)*lda];

      blasDAXPY( j-k, temp, work+k, 1, a+k+(j-1)*lda, 1 );
    }
  }

  return( info );
}

/**
 * Computes the lower-triangle Cholesky Decomposition
 * of a square, symmetric and positive-definite matrix.
 * Subroutines from LINPACK and BLAS are used.
 */
double **choleskyDecomposition( double **matrix, int n )
{
  int     i, j, k, info, *ipvt;
  double *a, *work, **result;
  
  a    = (double *) Malloc( n*n*sizeof( double ) );
  work = (double *) Malloc( n*sizeof( double ) );
  ipvt = (int *) Malloc( n*sizeof( int ) );

  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      a[k] = matrix[i][j];
      k++;
    }
    ipvt[i] = 0;
  }

  info = linpackDCHDC( a, n, n, work, ipvt );

  result = matrixNew( n, n );
  if( info != n ) /* Matrix is not positive definite */
  {
    k = 0;
    for( i = 0; i < n; i++ )
    {
      for( j = 0; j < n; j++ )
      {
        result[i][j] = i != j ? 0.0 : sqrt( matrix[i][j] );
        k++;
      }
    }
  }
  else
  {
    k = 0;
    for( i = 0; i < n; i++ )
    {
      for( j = 0; j < n; j++ )
      {
        result[i][j] = i < j ? 0.0 : a[k];
        k++;
      }
    }
  }

  free( ipvt );
  free( work );
  free( a );
  
  return( result );
}

/**
 * LINPACK subroutine.
 */
int linpackDTRDI( double t[], int ldt, int n )
{
  int    j, k, info;
  double temp;

  info = 0;
  for( k = n; 1 <= k; k-- )
  {
    if ( t[k-1+(k-1)*ldt] == 0.0 )
    {
      info = k;
      break;
    }

    t[k-1+(k-1)*ldt] = 1.0 / t[k-1+(k-1)*ldt];
    temp = -t[k-1+(k-1)*ldt];

    if ( k != n )
    {
      blasDSCAL( n-k, temp, t+k+(k-1)*ldt, 1 );
    }

    for( j = 1; j <= k-1; j++ )
    {
      temp = t[k-1+(j-1)*ldt];
      t[k-1+(j-1)*ldt] = 0.0;
      blasDAXPY( n-k+1, temp, t+k-1+(k-1)*ldt, 1, t+k-1+(j-1)*ldt, 1 );
    }
  }

  return( info );
}

/**
 * Computes the inverse of a matrix that is of
 * lower triangular form.
 */
double **matrixLowerTriangularInverse( double **matrix, int n )
{
  int     i, j, k, info;
  double *t, **result;
  
  t = (double *) Malloc( n*n*sizeof( double ) );

  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      t[k] = matrix[j][i];
      k++;
    }
  }

  info = linpackDTRDI( t, n, n );

  result = matrixNew( n, n );
  k = 0;
  for( i = 0; i < n; i++ )
  {
    for( j = 0; j < n; j++ )
    {
      result[j][i] = i > j ? 0.0 : t[k];
      k++;
    }
  }

  free( t );
  
  return( result );
}

/**
 * Writes the contents of a matrix of dimensions n0 x n1 to a file.
 */
void matrixWriteToFile( FILE *file, double **matrix, int n0, int n1 )
{
  int  i, j;
  char line_for_output[10000];
  
  sprintf( line_for_output, "[" );
  fputs( line_for_output, file );
  for( i = 0; i < n0; i++ )
  {
    sprintf( line_for_output, "[" );
    fputs( line_for_output, file );
    for( j = 0; j < n1; j++ )
    {
      sprintf( line_for_output, "%lf", matrix[i][j] );
      fputs( line_for_output, file );
      if( j < n1-1 )
      {
        sprintf( line_for_output, ", " );
        fputs( line_for_output, file );
      }
    }
    if( i == n0-1 )
      sprintf( line_for_output, "]" );
    else
      sprintf( line_for_output, "];" );
    fputs( line_for_output, file );
  }
  sprintf( line_for_output, "]\n" );
  fputs( line_for_output, file );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Merge Sort -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Sorts an array of doubles and returns the sort-order (small to large).
 */
int *mergeSort( double *array, int array_size )
{
  int i, *sorted, *tosort;

  sorted = (int *) Malloc( array_size * sizeof( int ) );
  tosort = (int *) Malloc( array_size * sizeof( int ) );
  for( i = 0; i < array_size; i++ )
    tosort[i] = i;

  if( array_size == 1 )
    sorted[0] = 0;
  else
    mergeSortWithinBounds( array, sorted, tosort, 0, array_size-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void mergeSortWithinBounds( double *array, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortWithinBounds( array, sorted, tosort, p, r );
    mergeSortWithinBounds( array, sorted, tosort, r+1, q );
    mergeSortMerge( array, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortMerge( double *array, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( array[tosort[i]] < array[tosort[j]] )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}

/**
 * Sorts an array of objectives and constraints
 * using constraint domination and returns the
 * sort-order (small to large).
 */
int *mergeSortFitness( double *objectives, double *constraints, int number_of_solutions )
{
  int i, *sorted, *tosort;

  sorted = (int *) Malloc( number_of_solutions * sizeof( int ) );
  tosort = (int *) Malloc( number_of_solutions * sizeof( int ) );
  for( i = 0; i < number_of_solutions; i++ )
    tosort[i] = i;

  if( number_of_solutions == 1 )
    sorted[0] = 0;
  else
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, 0, number_of_solutions-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the objectives and
 * constraints arrays between p and q.
 */
void mergeSortFitnessWithinBounds( double *objectives, double *constraints, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, p, r );
    mergeSortFitnessWithinBounds( objectives, constraints, sorted, tosort, r+1, q );
    mergeSortFitnessMerge( objectives, constraints, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortFitnessMerge( double *objectives, double *constraints, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( betterFitness( objectives[tosort[i]], constraints[tosort[i]],
                           objectives[tosort[j]], constraints[tosort[j]] ) )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=- Section Interpret Command Line -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Parses and checks the command line.
 */
void interpretCommandLine( int argc, char **argv )
{
  parseCommandLine( argc, argv );
  
  checkOptions();
}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void parseCommandLine( int argc, char **argv )
{
  int index;

  index = 1;

  parseOptions( argc, argv, &index );
  
  parseParameters( argc, argv, &index );
}

/**
 * Parses only the options from the command line.
 */
void parseOptions( int argc, char **argv, int *index )
{
  double dummy;

  write_generational_statistics = 0;
  write_generational_solutions  = 0;
  print_verbose_overview        = 0;
  use_vtr                       = 0;

  for( ; (*index) < argc; (*index)++ )
  {
    if( argv[*index][0] == '-' )
    {
      /* If it is a negative number, the option part is over */
      if( sscanf( argv[*index], "%lf", &dummy ) && argv[*index][1] != '\0' )
        break;

      if( argv[*index][1] == '\0' )
        optionError( argv, *index );
      else if( argv[*index][2] != '\0' )
        optionError( argv, *index );
      else
      {
        switch( argv[*index][1] )
        {
          case '?': printUsage(); break;
          case 'P': printAllInstalledProblems(); break;
          case 's': write_generational_statistics = 1; break;
          case 'w': write_generational_solutions  = 1; break;
          case 'v': print_verbose_overview        = 1; break;
          case 'r': use_vtr                       = 1; break;
          default : optionError( argv, *index );
        }
      }
    }
    else /* Argument is not an option, so option part is over */
     break;
  }
}

/**
 * Writes the names of all installed problems to the standard output.
 */
void printAllInstalledProblems( void )
{
  int i, n;

  n = numberOfInstalledProblems();
  printf("Installed optimization problems:\n");
  for( i = 0; i < n; i++ )
    printf("%3d: %s\n", i, installedProblemName( i ));

  exit( 0 );
}

/**
 * Informs the user of an illegal option and exits the program.
 */
void optionError( char **argv, int index )
{
  printf("Illegal option: %s\n\n", argv[index]);

  printUsage();
}

/**
 * Parses only the EA parameters from the command line.
 */
void parseParameters( int argc, char **argv, int *index )
{
  int noError;

  if( (argc - *index) != 9 )
  {
    printf("Number of parameters is incorrect, require 9 parameters (you provided %d).\n\n", (argc - *index));

    printUsage();
  }

  noError = 1;
  noError = noError && sscanf( argv[*index+0], "%d", &problem_index );
  noError = noError && sscanf( argv[*index+1], "%d", &number_of_parameters );
  noError = noError && sscanf( argv[*index+2], "%lf", &lower_user_range );
  noError = noError && sscanf( argv[*index+3], "%lf", &upper_user_range );
  noError = noError && sscanf( argv[*index+4], "%lf", &rotation_angle );
  noError = noError && sscanf( argv[*index+5], "%d", &maximum_number_of_evaluations );
  noError = noError && sscanf( argv[*index+6], "%lf", &vtr );
  noError = noError && sscanf( argv[*index+7], "%lf", &fitness_variance_tolerance );
  noError = noError && sscanf( argv[*index+8], "%d", &maximum_number_of_populations );

  if( !noError )
  {
    printf("Error parsing parameters.\n\n");

    printUsage();
  }
}

/**
 * Prints usage information and exits the program.
 */
void printUsage( void )
{
  printf("Usage: AMaLGaM-Full-Free [-?] [-P] [-s] [-w] [-v] [-r] pro dim low upp rot eva vtr tol pop\n");
  printf(" -?: Prints out this usage information.\n");
  printf(" -P: Prints out a list of all installed optimization problems.\n");
  printf(" -s: Enables computing and writing of statistics every generation.\n");
  printf(" -w: Enables writing of solutions and their fitnesses every generation.\n");
  printf(" -v: Enables verbose mode. Prints the settings before starting the run.\n");
  printf(" -r: Enables use of vtr in termination condition (value-to-reach).\n");
  printf("\n");
  printf("  pro: Index of optimization problem to be solved (minimization).\n");
  printf("  dim: Number of parameters.\n");
  printf("  low: Overall initialization lower bound.\n");
  printf("  upp: Overall initialization upper bound.\n");
  printf("  rot: The angle by which to rotate the problem.\n");
  printf("  eva: Maximum number of evaluations allowed.\n");
  printf("  vtr: The value to reach. If the objective value of the best feasible solution reaches\n");
  printf("       this value, termination is enforced (if -r is specified).\n");
  printf("  tol: The tolerance level for fitness variance (i.e. minimum fitness variance)\n");
  printf("  pop: The maximum number of parallel populations.\n");
  exit( 0 );
}

/**
 * Checks whether the selected options are feasible.
 */
void checkOptions( void )
{
  if( number_of_parameters < 1 )
  {
    printf("\n");
    printf("Error: number of parameters < 1 (read: %d). Require number of parameters >= 1.", number_of_parameters);
    printf("\n\n");

    exit( 0 );
  }

  if( maximum_number_of_populations < 1 )
  {
    printf("\n");
    printf("Error: maximum number of populations < 1 (read: %d). Require maximum number of populations >= 1.", maximum_number_of_populations);
    printf("\n\n");

    exit( 0 );
  }
  
  if( maximum_number_of_evaluations < 1 )
  {
    printf("\n");
    printf("Error: maximum number of evaluations < 1 (read: %d). Require maximum number of evaluations >= 1.", maximum_number_of_evaluations);
    printf("\n\n");

    exit( 0 );
  }
  
  if( installedProblemName( problem_index ) == NULL )
  {
    printf("\n");
    printf("Error: unknown index for problem (read index %d).", problem_index );
    printf("\n\n");

    exit( 0 );
  }
}

/**
 * Prints the settings as read on the command line.
 */
void printVerboseOverview( void )
{
  int i;
  
  printf("### Settings #################################################\n");
  printf("#\n");
  printf("# Statistics writing every generation: %s\n", write_generational_statistics ? "enabled" : "disabled");
  printf("# Population file writing            : %s\n", write_generational_solutions ? "enabled" : "disabled");
  printf("# Use of value-to-reach (vtr)        : %s\n", use_vtr ? "enabled" : "disabled");
  printf("#\n");
  printf("##############################################################\n");
  printf("#\n");
  printf("# Optimization problem    = %s\n", installedProblemName( problem_index ));
  printf("# Number of parameters    = %d\n", number_of_parameters);
  printf("# Initialization ranges   = ");
  for( i = 0; i < number_of_parameters; i++ )
  {
    printf("x_%d: [%e;%e]", i, lower_init_ranges[i], upper_init_ranges[i]);
    if( i < number_of_parameters-1 )
      printf("\n#                           ");
  }
  printf("\n");
  printf("# Boundary ranges         = ");
  for( i = 0; i < number_of_parameters; i++ )
  {
    printf("x_%d: [%e;%e]", i, lower_range_bounds[i], upper_range_bounds[i]);
    if( i < number_of_parameters-1 )
      printf("\n#                           ");
  }
  printf("\n");
  printf("# Rotation angle          = %e\n", rotation_angle);
  printf("# Maximum numb. of eval.  = %d\n", maximum_number_of_evaluations);
  printf("# Value to reach (vtr)    = %e\n", vtr);
  printf("# Fitness var. tolerance  = %e\n", fitness_variance_tolerance);
  printf("# Max. number of pop's    = %d\n", maximum_number_of_populations);
  printf("# Random seed             = %ld\n", random_seed);
  printf("#\n");
  printf("##############################################################\n");
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Random Numbers -=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns a random double, distributed uniformly between 0 and 1.
 */
double randomRealUniform01( void )
{
  int64_t n26, n27;
  double  result;

  random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n26                  = (int64_t)(random_seed_changing >> (48 - 26));
  random_seed_changing = (random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n27                  = (int64_t)(random_seed_changing >> (48 - 27));
  result               = (((int64_t)n26 << 27) + n27) / ((double) (1LLU << 53));

  return( result );
}
        
/**
 * Returns a random integer, distributed uniformly between 0 and maximum.
 */
int randomInt( int maximum )
{
  int result;
  
  result = (int) (((double) maximum)*randomRealUniform01());
  
  return( result );
}

/**
 * Returns a random double, distributed normally with mean 0 and variance 1.
 */
double random1DNormalUnit( void )
{
  double v1, v2, s, multiplier, value;

  if( haveNextNextGaussian )
  {
    haveNextNextGaussian = 0;

    return( nextNextGaussian );
  }
  else
  {
    do
    {
      v1 = 2 * (randomRealUniform01()) - 1;
      v2 = 2 * (randomRealUniform01()) - 1;
      s = v1 * v1 + v2 * v2;
    } while (s >= 1);

    value                = -2 * log(s)/s;
    multiplier           = value <= 0.0 ? 0.0 : sqrt( value );
    nextNextGaussian     = v2 * multiplier;
    haveNextNextGaussian = 1;

    return( v1 * multiplier );
  }
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Problems -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns the name of an installed problem.
 */
char *installedProblemName( int index )
{
  return( (char *) "Rosenbrock" );
}

/**
 * Returns the number of problems installed.
 */
int numberOfInstalledProblems( void )
{
  return 1;
}

/**
 * Returns the lower-range bound of an installed problem.
 */
double installedProblemLowerRangeBound( int index, int dimension )
{
  return( rosenbrockFunctionLowerRangeBound( dimension ) );
}

/**
 * Returns the upper-range bound of an installed problem.
 */
double installedProblemUpperRangeBound( int index, int dimension )
{
  return( rosenbrockFunctionUpperRangeBound( dimension ) );
}

/**
 * Returns whether a parameter is inside the range bound of
 * every problem.
 */
short isParameterInRangeBounds( double parameter, int dimension )
{
  if( parameter < installedProblemLowerRangeBound( problem_index, dimension ) ||
      parameter > installedProblemUpperRangeBound( problem_index, dimension ) ||
      isnan( parameter ) )
  {
    return( 0 );
  }
  
  return( 1 );
}

/**
 * Returns the value of the single objective
 * and the sum of all constraint violations
 * function after rotating the parameter vector.
 * Both are returned using pointer variables.
 */
void installedProblemEvaluation( int *pindex, double *parameters, double *objective_value, double *constraint_value, int *info )
{
  int index = *pindex;
  double *rotated_parameters;

  //number_of_evaluations++;

  if( rotation_angle == 0.0 )
    installedProblemEvaluationWithoutRotation( index, parameters, objective_value, constraint_value, info );
  else
  {
    rotated_parameters = matrixVectorMultiplication( rotation_matrix, parameters, number_of_parameters, number_of_parameters );

    installedProblemEvaluationWithoutRotation( index, rotated_parameters, objective_value, constraint_value, info );

    free( rotated_parameters );
  }
}

/**
 * Returns the value of the single objective
 * and the sum of all constraint violations
 * without rotating the parameter vector.
 * Both are returned using pointer variables.
 */
void installedProblemEvaluationWithoutRotation( int index, double *parameters, double *objective_value, double *constraint_value, int *info )
{
  *objective_value  = 0.0;
  *constraint_value = 0.0;

  rosenbrockFunctionProblemEvaluation( parameters, objective_value, constraint_value, info );
}

void rosenbrockFunctionProblemEvaluation( double *parameters, double *objective_value, double *constraint_value, int *info )
{
#if defined(RANDOMIZE)
  int    i;
  double result;
  
  result = 0.0;
  for( i = 0; i < number_of_parameters-1; i++ )
    result += 100*(parameters[i+1]-parameters[i]*parameters[i])*(parameters[i+1]-parameters[i]*parameters[i]) + (1.0-parameters[i])*(1.0-parameters[i]);

  usleep(1*1000);

  *objective_value  = result;
  *constraint_value = 0;
#else
  *objective_value  = -fitfun(parameters, number_of_parameters, NULL, info);
  *constraint_value = 0;
#endif
}

double rosenbrockFunctionLowerRangeBound( int dimension )
{
  return( -1e+308 );
}

double rosenbrockFunctionUpperRangeBound( int dimension )
{
  return( 1e+308 );
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initializations that are required before starting a run.
 */
void initialize( void )
{
  number_of_generations = 0;

  if( number_of_starts == 1 )
    number_of_evaluations = 0;

  alpha_AMS = 0.5*tau*(((double) population_size)/((double) (population_size-1)));
  delta_AMS = 2.0;

  initializeMemory();

  if( number_of_starts == 1 )
    initializeRandomNumberGenerator();

  initializeParameterRangeBounds();

  initializeDistributionMultipliers();

  initializeObjectiveRotationMatrix();

  initializePopulationsAndFitnessValues();

  computeRanks();
}

/**
 * Initializes the memory.
 */
void initializeMemory( void )
{
  int i, j;

  selection_size = (int) (tau*(population_size));

  populations                      = (double ***) Malloc( number_of_populations*sizeof( double ** ) );
  populations_terminated           = (short *) Malloc( number_of_populations*sizeof( short ) );
  no_improvement_stretch           = (int *) Malloc( number_of_populations*sizeof( int ) );
  objective_values                 = (double **) Malloc( number_of_populations*sizeof( double * ) );
  constraint_values                = (double **) Malloc( number_of_populations*sizeof( double * ) );
  ranks                            = (double **) Malloc( number_of_populations*sizeof( double * ) );
  selections                       = (double ***) Malloc( number_of_populations*sizeof( double ** ) );
  objective_values_selections      = (double **) Malloc( number_of_populations*sizeof( double * ) );
  constraint_values_selections     = (double **) Malloc( number_of_populations*sizeof( double * ) );
  mean_vectors                     = (double **) Malloc( number_of_populations*sizeof( double * ) );
  mean_vectors_previous            = (double **) Malloc( number_of_populations*sizeof( double * ) );
  covariance_matrices              = (double ***) Malloc( number_of_populations*sizeof( double ** ) );
  cholesky_factors_lower_triangle  = (double ***) Malloc( number_of_populations*sizeof( double ** ) );

  for( i = 0; i < number_of_populations; i++ )
  {
    populations[i] = (double **) Malloc( population_size*sizeof( double * ) );
    for( j = 0; j < population_size; j++ )
      populations[i][j] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    populations_terminated[i] = 0;

    no_improvement_stretch[i] = 0;

    objective_values[i] = (double *) Malloc( population_size*sizeof( double ) );
    
    constraint_values[i] = (double *) Malloc( population_size*sizeof( double ) );

    ranks[i] = (double *) Malloc( population_size*sizeof( double ) );

    selections[i] = (double **) Malloc( selection_size*sizeof( double * ) );
    for( j = 0; j < selection_size; j++ )
      selections[i][j] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    objective_values_selections[i] = (double *) Malloc( selection_size*sizeof( double ) );

    constraint_values_selections[i] = (double *) Malloc( selection_size*sizeof( double ) );

    mean_vectors[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    mean_vectors_previous[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    covariance_matrices[i] = (double **) Malloc( number_of_parameters*sizeof( double * ) );
    for( j = 0; j < number_of_parameters; j++ )
      covariance_matrices[i][j] = (double *) Malloc( number_of_parameters*sizeof( double ) );
    
    cholesky_factors_lower_triangle[i] = NULL;
  }

  lower_range_bounds = (double *) Malloc( number_of_parameters*sizeof( double ) );
  upper_range_bounds = (double *) Malloc( number_of_parameters*sizeof( double ) );
  lower_init_ranges  = (double *) Malloc( number_of_parameters*sizeof( double ) );
  upper_init_ranges  = (double *) Malloc( number_of_parameters*sizeof( double ) );
}


/**
 * Initializes the random number generator.
 */
void initializeRandomNumberGenerator( void )
{
#if 0
  struct tm *timep;
  time_t t;

  while( random_seed_changing == 0 )
  {
    t                    = time(NULL);
    timep                = localtime( &t );
    random_seed_changing = ((60*(long) timep->tm_min)) + (60*60*(long) timep->tm_hour) + ((long) timep->tm_sec);
    random_seed_changing = (random_seed_changing/((int) (9.99*randomRealUniform01())+1))*(((int) (randomRealUniform01()*1000000.0))%10);
  }

  random_seed = random_seed_changing;
#else
  random_seed = random_seed_changing = 280675;
#endif
}

/**
 * Initializes the parameter range bounds.
 */
void initializeParameterRangeBounds( void )
{
  int i;

  for( i = 0; i < number_of_parameters; i++ )
  {
    lower_range_bounds[i] = installedProblemLowerRangeBound( problem_index, i );
    upper_range_bounds[i] = installedProblemUpperRangeBound( problem_index, i );
  }

  for( i = 0; i < number_of_parameters; i++ )
  {
    lower_init_ranges[i] = lower_user_range;
    if( lower_user_range < lower_range_bounds[i] )
      lower_init_ranges[i] = lower_range_bounds[i];

    upper_init_ranges[i] = upper_user_range;
    if( upper_user_range > upper_range_bounds[i] )
      upper_init_ranges[i] = upper_range_bounds[i];
  }
}

/**
 * Initializes the distribution multipliers.
 */
void initializeDistributionMultipliers( void )
{
  int i;

  distribution_multipliers = (double *) Malloc( number_of_populations*sizeof( double ) );
  for( i = 0; i < number_of_populations; i++ )
    distribution_multipliers[i] = 1.0;

  samples_drawn_from_normal = (int *) Malloc( number_of_populations*sizeof( int ) );
  out_of_bounds_draws       = (int *) Malloc( number_of_populations*sizeof( int ) );

  distribution_multiplier_increase = 1.0/distribution_multiplier_decrease;
}

/**
 * Initializes the populations, the objective values and the constraint values.
 */
void initializePopulationsAndFitnessValues( void )
{
  int     i, j, k, o, q, *sorted, ssize, j_min, *temporary_population_sizes;
  double *distances, d, d_min, **solutions, *fitnesses, *constraints, **leader_vectors;

//  #pragma omp parallel
  {
//  #pragma omp single nowait
  {
  for( i = 0; i < number_of_populations; i++ )
  {
    for( j = 0; j < population_size; j++ )
    {
      for( k = 0; k < number_of_parameters; k++ )
        populations[i][j][k] = lower_init_ranges[k] + (upper_init_ranges[k] - lower_init_ranges[k])*randomRealUniform01();

//      #pragma omp task firstprivate(problem_index,i,j) shared(populations, objective_values, constraint_values)
      {
	int info[5];
	info[0] = number_of_starts;		// start
	info[1] = 0;				// initialize
	info[2] = number_of_generations;	// generation (step)
	info[3] = i;				// population
	info[4] = j;				// element
//      installedProblemEvaluation( problem_index, populations[i][j], &(objective_values[i][j]), &(constraint_values[i][j]) );
	torc_create(-1, installedProblemEvaluation, 5,
		1, MPI_INT, CALL_BY_COP,
		number_of_parameters, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_RES,
		1, MPI_DOUBLE, CALL_BY_RES,
		5, MPI_INT, CALL_BY_COP,
		 &problem_index, populations[i][j], &(objective_values[i][j]), &(constraint_values[i][j]), info);
	number_of_evaluations++;
      }
    }
  }
//  #pragma omp taskwait
	torc_waitall();
  }
  }

  /* Initialize means and redistribute solutions */
  if( number_of_populations > 1 )
  {
    ssize                      = number_of_populations*population_size;
    solutions                  = (double **) Malloc( ssize*sizeof( double * ) );
    fitnesses                  = (double *) Malloc( ssize*sizeof( double ) );
    constraints                = (double *) Malloc( ssize*sizeof( double ) );
    temporary_population_sizes = (int *) Malloc( ssize*sizeof( int ) );
    distances                  = (double *) Malloc( ssize*sizeof( double ) );
    leader_vectors             = (double **) Malloc( number_of_populations*sizeof( double * ) );

    for( i = 0; i < ssize; i++ )
      solutions[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    for( i = 0; i < number_of_populations; i++ )
      leader_vectors[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

    q = 0;
    for( i = 0; i < number_of_populations; i++ )
    {
      for( j = 0; j < population_size; j++ )
      {
        for( k = 0; k < number_of_parameters; k++ )
          solutions[q][k] = populations[i][j][k];

        fitnesses[q]   = objective_values[i][j];
        constraints[q] = constraint_values[i][j];
        q++;
      }
    }
    o = randomInt( number_of_parameters );
  
    for( i = 0; i < ssize; i++ )
      distances[i] = solutions[i][o];
  
    sorted = mergeSort( distances, ssize );
  
    q = 0;
    for( i = 0; i < number_of_parameters; i++ )
      leader_vectors[q][i] = solutions[sorted[0]][i];
  
    for( i = 0; i < ssize; i++ )
      distances[i] = -distanceInParameterSpace( leader_vectors[q], solutions[i] );

    free( sorted );
  
    q++;
    while( q < number_of_populations )
    {
      sorted = mergeSort( distances, ssize );
  
      for( i = 0; i < number_of_parameters; i++ )
        leader_vectors[q][i] = solutions[sorted[0]][i];
  
      for( i = 0; i < ssize; i++ )
      {
        d = -distanceInParameterSpace( leader_vectors[q], solutions[i] );
        if( d > distances[i] )
          distances[i] = d;
      }
  
      free( sorted );

      q++;
    }
  
    for( i = 0; i < number_of_populations; i++ )
      temporary_population_sizes[i] = 0;
  
    for( i = 0; i < ssize; i++ )
    {
      j_min = -1;
      d_min = 0;
      for( j = 0; j < number_of_populations; j++ )
      {
        if( temporary_population_sizes[j] < population_size )
        {
          d = distanceInParameterSpace( solutions[i], leader_vectors[j] );
          if( (j_min == -1) || (d < d_min) )
          {
            j_min = j;
            d_min = d;
          }
        }
      }
      for( k = 0; k < number_of_parameters; k++ )
        populations[j_min][temporary_population_sizes[j_min]][k]  = solutions[i][k];
      objective_values[j_min][temporary_population_sizes[j_min]]  = fitnesses[i];
      constraint_values[j_min][temporary_population_sizes[j_min]] = constraints[i];
      temporary_population_sizes[j_min]++;
    }

    for( i = 0; i < number_of_populations; i++ )
      free( leader_vectors[i] );
    free( leader_vectors );
    free( distances );
    free( temporary_population_sizes );
    free( fitnesses );
    free( constraints );
    for( i = 0; i < ssize; i++ )
      free( solutions[i] );
    free( solutions );
  }
}

/**
 * Computes the rotation matrix to be applied to any solution
 * before evaluating it (i.e. turns the evaluation functions
 * into rotated evaluation functions.
 */
void initializeObjectiveRotationMatrix( void )
{
  int      i, j, index0, index1;
  double **matrix, **product, theta, cos_theta, sin_theta;

  if( rotation_angle == 0.0 )
    return;

  matrix = (double **) Malloc( number_of_parameters*sizeof( double * ) );
  for( i = 0; i < number_of_parameters; i++ )
    matrix[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

  rotation_matrix = (double **) Malloc( number_of_parameters*sizeof( double * ) );
  for( i = 0; i < number_of_parameters; i++ )
    rotation_matrix[i] = (double *) Malloc( number_of_parameters*sizeof( double ) );

  /* Initialize the rotation matrix to the identity matrix */
  for( i = 0; i < number_of_parameters; i++ )
  {
    for( j = 0; j < number_of_parameters; j++ )
      rotation_matrix[i][j] = 0.0;
    rotation_matrix[i][i] = 1.0;
  }

  /* Construct all rotation matrices (quadratic number) and multiply */
  theta     = (rotation_angle/180.0)*PI;
  cos_theta = cos( theta );
  sin_theta = sin( theta );
  for( index0 = 0; index0 < number_of_parameters-1; index0++ )
  {
    for( index1 = index0+1; index1 < number_of_parameters; index1++ )
    {
      for( i = 0; i < number_of_parameters; i++ )
      {
        for( j = 0; j < number_of_parameters; j++ )
          matrix[i][j] = 0.0;
        matrix[i][i] = 1.0;
      }
      matrix[index0][index0] = cos_theta;
      matrix[index0][index1] = -sin_theta;
      matrix[index1][index0] = sin_theta;
      matrix[index1][index1] = cos_theta;
      
      product = matrixMatrixMultiplication( matrix, rotation_matrix, number_of_parameters, number_of_parameters, number_of_parameters );
      for( i = 0; i < number_of_parameters; i++ )
        for( j = 0; j < number_of_parameters; j++ )
          rotation_matrix[i][j] = product[i][j];

      for( i = 0; i < number_of_parameters; i++ )
        free( product[i] );
      free( product );
    }
  }

  for( i = 0; i < number_of_parameters; i++ )
    free( matrix[i] );
  free( matrix );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Ranking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Computes the ranks of all populations.
 */
void computeRanks( void )
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
    if( !populations_terminated[i] )
      computeRanksForOnePopulation( i );
}

/**
 * Computes the ranks for one population.
 */
void computeRanksForOnePopulation( int population_index )
{
  int i, *sorted, rank;

  sorted = mergeSortFitness( objective_values[population_index], constraint_values[population_index], population_size );
  
  rank                               = 0;
  ranks[population_index][sorted[0]] = rank;
  for( i = 1; i < population_size; i++ )
  {
    if( objective_values[population_index][sorted[i]] != objective_values[population_index][sorted[i-1]] )
      rank++;

    ranks[population_index][sorted[i]] = rank;
  }

  free( sorted );
}

/**
 * Computes the distance between two solutions a and b as
 * the Euclidean distance in parameter space.
 */
double distanceInParameterSpace( double *solution_a, double *solution_b )
{
  int    i;
  double value, result;
  
  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    value   = solution_b[i] - solution_a[i];
    result += value*value;
  }
  result = sqrt( result );
  
  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Output =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Writes (appends) statistics about the current generation to a
 * file named "statistics.dat".
 */
void writeGenerationalStatistics( void )
{
  int     i, j;
  char    string[1000];
  double  overall_objective_avg, overall_objective_var, overall_objective_best, overall_objective_worst,
          overall_constraint_avg, overall_constraint_var, overall_constraint_best, overall_constraint_worst,
         *population_objective_avg, *population_objective_var, *population_objective_best, *population_objective_worst,
         *population_constraint_avg, *population_constraint_var, *population_constraint_best, *population_constraint_worst;
  FILE   *file;

  /* First compute the statistics */
  population_objective_avg    = (double *) Malloc( number_of_populations*sizeof( double ) );
  population_constraint_avg   = (double *) Malloc( number_of_populations*sizeof( double ) );
  population_objective_var    = (double *) Malloc( number_of_populations*sizeof( double ) );
  population_constraint_var   = (double *) Malloc( number_of_populations*sizeof( double ) );
  population_objective_best   = (double *) Malloc( number_of_populations*sizeof( double ) );
  population_constraint_best  = (double *) Malloc( number_of_populations*sizeof( double ) );
  population_objective_worst  = (double *) Malloc( number_of_populations*sizeof( double ) );
  population_constraint_worst = (double *) Malloc( number_of_populations*sizeof( double ) );

  /* Overall */
  /* Average, best and worst */
  overall_objective_avg    = 0.0;
  overall_constraint_avg   = 0.0;
  overall_objective_best   = objective_values[0][0];
  overall_objective_worst  = objective_values[0][0];
  overall_constraint_best  = constraint_values[0][0];
  overall_constraint_worst = constraint_values[0][0];
  for( i = 0; i < number_of_populations; i++ )
  {
    for( j = 0; j < population_size; j++ )
    {
      overall_objective_avg += objective_values[i][j];
      overall_constraint_avg += constraint_values[i][j];
      if( betterFitness( objective_values[i][j], constraint_values[i][j], overall_objective_best, overall_constraint_best ) )
      {
        overall_objective_best  = objective_values[i][j];
        overall_constraint_best = constraint_values[i][j];
      }
      if( betterFitness( overall_objective_worst, overall_constraint_worst, objective_values[i][j], constraint_values[i][j] ) )
      {
        overall_objective_worst  = objective_values[i][j];
        overall_constraint_worst = constraint_values[i][j];
      }
    }
  }
  overall_objective_avg = overall_objective_avg / ((double) (number_of_populations*population_size));
  overall_constraint_avg = overall_constraint_avg / ((double) (number_of_populations*population_size));

  /* Variance */
  overall_objective_var    = 0.0;
  overall_constraint_var   = 0.0;
  for( i = 0; i < number_of_populations; i++ )
  {
    for( j = 0; j < population_size; j++ )
    {
      overall_objective_var += (objective_values[i][j] - overall_objective_avg)*(objective_values[i][j] - overall_objective_avg);
      overall_constraint_var += (constraint_values[i][j] - overall_constraint_avg)*(constraint_values[i][j] - overall_constraint_avg);
    }
  }
  overall_objective_var = overall_objective_var / ((double) (number_of_populations*population_size));
  overall_constraint_var = overall_constraint_var / ((double) (number_of_populations*population_size));

  if( overall_objective_var <= 0.0 )
     overall_objective_var = 0.0;
  if( overall_constraint_var <= 0.0 )
     overall_constraint_var = 0.0;

  /* Per population */
  for( i = 0; i < number_of_populations; i++ )
  {
    /* Average, best and worst */
    population_objective_avg[i]    = 0.0;
    population_constraint_avg[i]   = 0.0;
    population_objective_best[i]   = objective_values[i][0];
    population_constraint_best[i]  = constraint_values[i][0];
    population_objective_worst[i]  = objective_values[i][0];
    population_constraint_worst[i] = constraint_values[i][0];
    for( j = 0; j < population_size; j++ )
    {
      population_objective_avg[i]  += objective_values[i][j];
      population_constraint_avg[i] += constraint_values[i][j];
      if( betterFitness( objective_values[i][j], constraint_values[i][j], population_objective_best[i], population_constraint_best[i] ) )
      {
        population_objective_best[i] = objective_values[i][j];
        population_constraint_best[i] = constraint_values[i][j];
      }
      if( betterFitness( population_objective_worst[i], population_constraint_worst[i], objective_values[i][j], constraint_values[i][j] ) )
      {
        population_objective_worst[i] = objective_values[i][j];
        population_constraint_worst[i] = constraint_values[i][j];
      }
    }
    population_objective_avg[i]  = population_objective_avg[i] / ((double) population_size);
    population_constraint_avg[i] = population_constraint_avg[i] / ((double) population_size);

    /* Variance */
    population_objective_var[i]    = 0.0;
    population_constraint_var[i]   = 0.0;
    for( j = 0; j < population_size; j++ )
    {
      population_objective_var[i]  += (objective_values[i][j] - population_objective_avg[i])*(objective_values[i][j] - population_objective_avg[i]);
      population_constraint_var[i] += (constraint_values[i][j] - population_constraint_avg[i])*(constraint_values[i][j] - population_constraint_avg[i]);
    }
    population_objective_var[i]  = population_objective_var[i] / ((double) population_size);
    population_constraint_var[i] = population_constraint_var[i] / ((double) population_size);

    if( population_objective_var[i] <= 0.0 )
       population_objective_var[i] = 0.0;
    if( population_constraint_var[i] <= 0.0 )
       population_constraint_var[i] = 0.0;
  }

  /* Then write them */
  file = NULL;
  if( number_of_generations == 0 && number_of_starts == 1 )
  {
    file = fopen( "statistics.dat", "w" );

    sprintf( string, "# Generation Evaluations  Best-obj-ever Best-con-ever Average-obj. Variance-obj.     Best-obj.    Worst-obj.  Average-con. Variance-con.     Best-con.    Worst-con.   [ ");
    fputs( string, file );

    for( i = 0; i < number_of_populations; i++ )
    {
      sprintf( string, "Pop.index     Dis.mult.  Pop.avg.obj.  Pop.var.obj. Pop.best.obj. Pop.worst.obj.  Pop.avg.con.  Pop.var.con. Pop.best.con. Pop.worst.con." );
      fputs( string, file );
      if( i < number_of_populations-1 )
      {
        sprintf( string, " | " );
        fputs( string, file );
      }
    }
    sprintf( string, " ]\n" );
    fputs( string, file );
  }
  else
    file = fopen( "statistics.dat", "a" );

  if( number_of_starts == 1 ||
      betterFitness( overall_objective_best,
                     overall_constraint_best,
                     best_so_far_objective_value,
                     best_so_far_constraint_value ) )
    sprintf( string, "  %10d %11d %13e %13e %13e %13e %13e %13e %13e %13e %13e %13e   [ ", number_of_generations, number_of_evaluations, overall_objective_best, overall_constraint_best, overall_objective_avg, overall_objective_var, overall_objective_best, overall_objective_worst, overall_constraint_avg, overall_constraint_var, overall_constraint_best, overall_constraint_worst );
  else
    sprintf( string, "  %10d %11d %13e %13e %13e %13e %13e %13e %13e %13e %13e %13e   [ ", number_of_generations, number_of_evaluations, best_so_far_objective_value, best_so_far_constraint_value, overall_objective_avg, overall_objective_var, overall_objective_best, overall_objective_worst, overall_constraint_avg, overall_constraint_var, overall_constraint_best, overall_constraint_worst );
  fputs( string, file );

  for( i = 0; i < number_of_populations; i++ )
  {
    sprintf( string, "%9d %13e %13e %13e %13e  %13e %13e %13e %13e  %13e", i, distribution_multipliers[i], population_objective_avg[i], population_objective_var[i], population_objective_best[i], population_objective_worst[i], population_constraint_avg[i], population_constraint_var[i], population_constraint_best[i], population_constraint_worst[i] );
    fputs( string, file );
    if( i < number_of_populations-1 )
    {
      sprintf( string, " | " );
      fputs( string, file );
    }
  }
  sprintf( string, " ]\n" );
  fputs( string, file );
  
  fclose( file );

  free( population_objective_avg );
  free( population_constraint_avg );
  free( population_objective_var );
  free( population_constraint_var );
  free( population_objective_best );
  free( population_constraint_best );
  free( population_objective_worst );
  free( population_constraint_worst );
}

/**
 * Writes the solutions to various files. The filenames
 * contain the generation. If the flag final is set
 * (final != 0), the generation number in the filename
 * is replaced with the word "final".
 *
 * all_populations_generation_xxxxx.dat : all populations combined
 * population_xxxxx_generation_xxxxx.dat: the individual populations
 * selection_xxxxx_generation_xxxxx.dat : the individual selections
 */
void writeGenerationalSolutions( short final )
{
  int   i, j, k;
  char  string[1000];
  FILE *file_all, *file_population, *file_selection;

  if( final )
    sprintf( string, "all_populations_generation_final.dat" );
  else
    sprintf( string, "all_populations_generation_%05d.dat", number_of_generations );
  file_all = fopen( string, "w" );
  
  for( i = 0; i < number_of_populations; i++ )
  {
    if( final )
      sprintf( string, "population_%05d_generation_final.dat", i );
    else
      sprintf( string, "population_%05d_generation_%05d.dat", i, number_of_generations );
    file_population = fopen( string, "w" );

    if( number_of_generations > 0 && !final )
    {
      sprintf( string, "selection_%05d_generation_%05d.dat", i, number_of_generations-1 );
      file_selection = fopen( string, "w" );
    }
    
    /* Populations */
    for( j = 0; j < population_size; j++ )
    {
      for( k = 0; k < number_of_parameters; k++ )
      {
        sprintf( string, "%13e", populations[i][j][k] );
        fputs( string, file_all );
        fputs( string, file_population );
        if( k < number_of_parameters-1 )
        {
          sprintf( string, " " );
          fputs( string, file_all );
          fputs( string, file_population );
        }
      }
      sprintf( string, "     " );
      fputs( string, file_all );
      fputs( string, file_population );
      sprintf( string, "%13e %13e", objective_values[i][j], constraint_values[i][j] );
      fputs( string, file_all );
      fputs( string, file_population );
      sprintf( string, "\n" );
      fputs( string, file_all );
      fputs( string, file_population );
    }
    
    fclose( file_population );
    
    /* Selections */
    if( number_of_generations > 0 && !final )
    {
      for( j = 0; j < selection_size; j++ )
      {
        for( k = 0; k < number_of_parameters; k++ )
        {
         sprintf( string, "%13e", selections[i][j][k] );
         fputs( string, file_selection );
         if( k < number_of_parameters-1 )
         {
           sprintf( string, " " );
           fputs( string, file_selection );
         }
         sprintf( string, "     " );
         fputs( string, file_selection );
        }
        sprintf( string, "%13e %13e", objective_values_selections[i][j], constraint_values_selections[i][j] );
        fputs( string, file_selection );
        sprintf( string, "\n" );
        fputs( string, file_selection );
      }
      fclose( file_selection );
    }
  }
  
  fclose( file_all );

  writeGenerationalSolutionsBest( final );
}

/**
 * Writes the best solution (measured in the single
 * available objective) to a file named
 * best_generation_xxxxx.dat where xxxxx is the
 * generation number. If the flag final is set
 * (final != 0), the generation number in the filename
 * is replaced with the word "final". The output
 * file contains the solution values with the
 * dimensions separated by a single white space,
 * followed by five white spaces and then the
 * single objective value for that solution
 * and its sum of constraint violations.
 */
void writeGenerationalSolutionsBest( short final )
{
  if( final )
    writeGenerationalSolutionsBestFinal();
  else
    writeGenerationalSolutionsBestNotFinal();
}

/**
 * Writes the best solution upon termination.
 * The best solution ever found is used to this end.
 */
void writeGenerationalSolutionsBestFinal( void )
{
  int   i;
  char  string[1000];
  FILE *file;

  sprintf( string, "best_generation_final.dat" );
  file = fopen( string, "w" );

  for( i = 0; i < number_of_parameters; i++ )
  {
    sprintf( string, "%13e", best_so_far_solution[i] );
    fputs( string, file );
    if( i < number_of_parameters-1 )
    {
      sprintf( string, " " );
      fputs( string, file );
    }
  }
  sprintf( string, "     " );
  fputs( string, file );
  sprintf( string, "%13e %13e", best_so_far_objective_value, best_so_far_constraint_value );
  fputs( string, file );
  sprintf( string, "\n" );
  fputs( string, file );

  fclose( file );
}

/**
 * Writes the best solution while the algorithm
 * is running. The best solution in all populations
 * is determined to this end.
 */
void writeGenerationalSolutionsBestNotFinal( void )
{
  int   i, population_index_best, individual_index_best;
  char  string[1000];
  FILE *file;

  /* First find the best of all */
  determineBestSolutionInCurrentPopulations( &population_index_best, &individual_index_best );
  
  /* Then output it */
  sprintf( string, "best_generation_%05d.dat", number_of_generations );
  file = fopen( string, "w" );

  for( i = 0; i < number_of_parameters; i++ )
  {
    sprintf( string, "%13e", populations[population_index_best][individual_index_best][i] );
    fputs( string, file );
    if( i < number_of_parameters-1 )
    {
      sprintf( string, " " );
      fputs( string, file );
    }
  }
  sprintf( string, "     " );
  fputs( string, file );
  sprintf( string, "%13e %13e", objective_values[population_index_best][individual_index_best], constraint_values[population_index_best][individual_index_best] );
  fputs( string, file );
  sprintf( string, "\n" );
  fputs( string, file );

  fclose( file );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Termination -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns 1 if termination should be enforced
 * for the overall run with increasing population
 * size and increasing number of parallel populations,
 * 0 otherwise.
 */
short checkTerminationCondition( void )
{
  if( checkNumberOfEvaluationsTerminationCondition() )
    return( 1 );
  
  if( use_vtr )
  {
    if( best_so_far_constraint_value == 0 && best_so_far_objective_value <= vtr  )
      return( 1 );
    
    return( 0 );
  }

  return( 0 );
}

/**
 * Returns 1 if termination should be enforced
 * for a single run, 0 otherwise.
 */
short checkTerminationConditionForRunOnce( void )
{
  short allTrue;
  int   i;

  if( checkNumberOfEvaluationsTerminationCondition() )
    return( 1 );
  
  if( use_vtr )
  {
    if( checkVTRTerminationCondition() )
      return( 1 );
  }

  checkFitnessVarianceTermination();

  checkDistributionMultiplierTerminationCondition();

  allTrue = 1;
  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      allTrue = 0;
      break;
    }
  }

  return( allTrue );
}

/**
 * Returns 1 if the maximum number of evaluations
 * has been reached, 0 otherwise.
 */
short checkNumberOfEvaluationsTerminationCondition( void )
{
  if( number_of_evaluations >= maximum_number_of_evaluations )
    return( 1 );

  return( 0 );
}

/**
 * Returns 1 if the value-to-reach has been reached (in any population).
 */
short checkVTRTerminationCondition( void )
{
  int population_of_best, index_of_best;

  determineBestSolutionInCurrentPopulations( &population_of_best, &index_of_best );

  if( constraint_values[population_of_best][index_of_best] == 0 && objective_values[population_of_best][index_of_best] <= vtr  )
    return( 1 );
    
  return( 0 );
}

/**
 * Determines which solution is the best of all solutions
 * in all current populations.
 */
void determineBestSolutionInCurrentPopulations( int *population_of_best, int *index_of_best )
{
  int i, j;

  (*population_of_best) = 0;
  (*index_of_best)      = 0;
  for( i = 0; i < number_of_populations; i++ )
  {
    for( j = 0; j < population_size; j++ )
    {
      if( betterFitness( objective_values[i][j], constraint_values[i][j],
                         objective_values[(*population_of_best)][(*index_of_best)], constraint_values[(*population_of_best)][(*index_of_best)] ) )
      {
        (*population_of_best) = i;
        (*index_of_best)      = j;
      }
    }
  }
}

void checkFitnessVarianceTermination( void )
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
      if( checkFitnessVarianceTerminationSinglePopulation( i ) )
        populations_terminated[i] = 1;
  }
}

/**
 * Returns 1 if the fitness variance in a specific population
 * has become too small (user-defined tolerance).
 */
short checkFitnessVarianceTerminationSinglePopulation( int population_index )
{
  int    i;
  double objective_avg, objective_var;
  
  objective_avg = 0.0;
  for( i = 0; i < population_size; i++ )
    objective_avg  += objective_values[population_index][i];
  objective_avg = objective_avg / ((double) population_size);

  objective_var = 0.0;
  for( i = 0; i < population_size; i++ )
    objective_var  += (objective_values[population_index][i]-objective_avg)*(objective_values[population_index][i]-objective_avg);
  objective_var = objective_var / ((double) population_size);

  if( objective_var <= 0.0 )
    objective_var = 0.0;

  if( objective_var <= fitness_variance_tolerance )
    return( 1 );

  return( 0 );
}

/**
 * Checks whether the distribution multiplier in any population
 * has become too small (1e-10).
 */
void checkDistributionMultiplierTerminationCondition( void )
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
      if( distribution_multipliers[i] < 1e-10 )
        populations_terminated[i] = 1;
  }
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Selection =-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Makes a set of selected solutions for each population.
 */
void makeSelections( void )
{
  int i;
  
  for( i = 0; i < number_of_populations; i++ )
    if( !populations_terminated[i] )
      makeSelectionsForOnePopulation( i );
}

/**
 * Performs truncation selection on a single population.
 */
void makeSelectionsForOnePopulation( int population_index )
{
  int i, j, *sorted;
  
  sorted = mergeSort( ranks[population_index], population_size );

  if( ranks[population_index][sorted[selection_size-1]] == 0 )
    makeSelectionsForOnePopulationUsingDiversityOnRank0( population_index );
  else
  {
    for( i = 0; i < selection_size; i++ )
    {
      for( j = 0; j < number_of_parameters; j++ )
        selections[population_index][i][j] = populations[population_index][sorted[i]][j];

      objective_values_selections[population_index][i]  = objective_values[population_index][sorted[i]];
      constraint_values_selections[population_index][i] = constraint_values[population_index][sorted[i]];
    }
  }
  
  free( sorted );
}

/**
 * Performs selection from all solutions that have rank 0
 * based on diversity.
 */
void makeSelectionsForOnePopulationUsingDiversityOnRank0( int population_index )
{
  int     i, j, number_of_rank0_solutions, *preselection_indices,
         *selection_indices, index_of_farthest, number_selected_so_far;
  double *nn_distances, distance_of_farthest, value;

  number_of_rank0_solutions = 0;
  for( i = 0; i < population_size; i++ )
  {
    if( ranks[population_index][i] == 0 )
      number_of_rank0_solutions++;
  }

  preselection_indices = (int *) Malloc( number_of_rank0_solutions*sizeof( int ) );
  j                    = 0;
  for( i = 0; i < population_size; i++ )
  {
    if( ranks[population_index][i] == 0 )
    {
      preselection_indices[j] = i;
      j++;
    }
  }

  index_of_farthest    = 0;
  distance_of_farthest = objective_values[population_index][preselection_indices[0]];
  for( i = 1; i < number_of_rank0_solutions; i++ )
  {
    if( objective_values[population_index][preselection_indices[i]] > distance_of_farthest )
    {
      index_of_farthest    = i;
      distance_of_farthest = objective_values[population_index][preselection_indices[i]];
    }
  }

  number_selected_so_far                    = 0;
  selection_indices                         = (int *) Malloc( selection_size*sizeof( int ) );
  selection_indices[number_selected_so_far] = preselection_indices[index_of_farthest];
  preselection_indices[index_of_farthest]   = preselection_indices[number_of_rank0_solutions-1];
  number_of_rank0_solutions--;
  number_selected_so_far++;

  nn_distances = (double *) Malloc( number_of_rank0_solutions*sizeof( double ) );
  for( i = 0; i < number_of_rank0_solutions; i++ )
    nn_distances[i] = distanceInParameterSpace( populations[population_index][preselection_indices[i]], populations[population_index][selection_indices[number_selected_so_far-1]] );
  
  while( number_selected_so_far < selection_size )
  {
    index_of_farthest    = 0;
    distance_of_farthest = nn_distances[0];
    for( i = 1; i < number_of_rank0_solutions; i++ )
    {
      if( nn_distances[i] > distance_of_farthest )
      {
        index_of_farthest    = i;
        distance_of_farthest = nn_distances[i];
      }
    }
    
    selection_indices[number_selected_so_far] = preselection_indices[index_of_farthest];
    preselection_indices[index_of_farthest]   = preselection_indices[number_of_rank0_solutions-1];
    nn_distances[index_of_farthest]           = nn_distances[number_of_rank0_solutions-1];
    number_of_rank0_solutions--;
    number_selected_so_far++;

    for( i = 0; i < number_of_rank0_solutions; i++ )
    {
      value = distanceInParameterSpace( populations[population_index][preselection_indices[i]], populations[population_index][selection_indices[number_selected_so_far-1]] );
      if( value < nn_distances[i] )
        nn_distances[i] = value;
    }
  }

  for( i = 0; i < selection_size; i++ )
  {
    for( j = 0; j < number_of_parameters; j++ )
      selections[population_index][i][j] = populations[population_index][selection_indices[i]][j];

    objective_values_selections[population_index][i]  = objective_values[population_index][selection_indices[i]];
    constraint_values_selections[population_index][i] = constraint_values[population_index][selection_indices[i]];
  }

  free( nn_distances );
  free( selection_indices );
  free( preselection_indices );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Variation -==-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * First estimates the parameters of a normal distribution in the
 * parameter space from the selected sets of solutions (a separate
 * normal distribution for each population). Then copies the single
 * best selected solutions to their respective populations. Finally
 * fills up each population, after the variances have been scaled,
 * by drawing new samples from the normal distributions and applying
 * AMS to several of these new solutions. Then, the fitness ranks
 * are recomputed. Finally, the distribution multipliers are adapted
 * according to the SDR-AVS mechanism.
 */
void makePopulations( void )
{
  estimateParametersAllPopulations();
  
  copyBestSolutionsToPopulations();

  applyDistributionMultipliers();

  generateAndEvaluateNewSolutionsToFillPopulations();

  computeRanks();

  adaptDistributionMultipliers();
}

/**
 * Estimates the parameters of the multivariate normal
 * distribution for each population separately.
 */
void estimateParametersAllPopulations( void )
{
  int i;

  for( i = 0; i < number_of_populations; i++ )
    if( !populations_terminated[i] )
      estimateParameters( i );
}

/**
 * Estimates the paramaters of the multivariate
 * normal distribution for a specified population.
 */
void estimateParameters( int population_index )
{
  estimateParametersML( population_index );
}

/**
 * Estimates (with maximum likelihood) the
 * parameters of a multivariate normal distribution
 * for a specified population.
 */
void estimateParametersML( int population_index )
{
  estimateMeanVectorML( population_index );

  estimateCovarianceMatrixML( population_index );
}

/**
 * Computes the sample mean for a specified population.
 */
void estimateMeanVectorML( int population_index )
{
  int i, j;

  if( number_of_generations > 0 )
  {
    for( i = 0; i < number_of_parameters; i++ )
      mean_vectors_previous[population_index][i] = mean_vectors[population_index][i];
  }

  for( i = 0; i < number_of_parameters; i++ )
  {
    mean_vectors[population_index][i] = 0.0;

    for( j = 0; j < selection_size; j++ )
      mean_vectors[population_index][i] += selections[population_index][j][i];

    mean_vectors[population_index][i] /= (double) selection_size;
  }

  /* Change the focus of the search to the best solution */
  if( distribution_multipliers[population_index] < 1.0 )
    for( i = 0; i < number_of_parameters; i++ )
      mean_vectors[population_index][i] = selections[population_index][0][i];
}

/**
 * Computes the matrix of sample covariances for
 * a specified population.
 *
 * It is important that the pre-condition must be satisified:
 * estimateMeanVector was called first.
 */
void estimateCovarianceMatrixML( int population_index )
{
  int i, j, k;

  /* First do the maximum-likelihood estimate from data */
  for( i = 0; i < number_of_parameters; i++ )
  {
    for( j = i; j < number_of_parameters; j++ )
    {
      covariance_matrices[population_index][i][j] = 0.0;

      for( k = 0; k < selection_size; k++ )
        covariance_matrices[population_index][i][j] += (selections[population_index][k][i]-mean_vectors[population_index][i])*(selections[population_index][k][j]-mean_vectors[population_index][j]);

      covariance_matrices[population_index][i][j] /= (double) selection_size;
    }
  }

  for( i = 0; i < number_of_parameters; i++ )
    for( j = 0; j < i; j++ )
      covariance_matrices[population_index][i][j] = covariance_matrices[population_index][j][i];
}

/**
 * Copies the single very best of the selected solutions
 * to their respective populations.
 */
void copyBestSolutionsToPopulations( void )
{
  int i, k;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      for( k = 0; k < number_of_parameters; k++ )
        populations[i][0][k] = selections[i][0][k];

      objective_values[i][0]  = objective_values_selections[i][0];
      constraint_values[i][0] = constraint_values_selections[i][0];
    }
  }
}

/**
 * Applies the distribution multipliers.
 */
void applyDistributionMultipliers( void )
{
  int i, j, k;
  
  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      for( j = 0; j < number_of_parameters; j++ )
        for( k = 0; k < number_of_parameters; k++ )
          covariance_matrices[i][j][k] *= distribution_multipliers[i];
    }
  }
}

/**
 * Generates new solutions for each
 * of the populations in turn.
 */
void generateAndEvaluateNewSolutionsToFillPopulations( void )
{
  short   out_of_range;
  int     i, j, k, q, number_of_AMS_solutions;
  double *solution, *solution_AMS, shrink_factor;

  solution_AMS = (double *) Malloc( number_of_parameters*sizeof( double ) );

//  #pragma omp parallel
  {
//  #pragma omp single nowait
  {
  for( i = 0; i < number_of_populations; i++ )
  {
    computeParametersForSampling( i );

    if( !populations_terminated[i] )
    {
      number_of_AMS_solutions      = (int) (alpha_AMS*(population_size-1));
      samples_drawn_from_normal[i] = 0;
      out_of_bounds_draws[i]       = 0;
      q                            = 0;
      for( j = 1; j < population_size; j++ )
      {
        solution = generateNewSolution( i );
  
        for( k = 0; k < number_of_parameters; k++ )
          populations[i][j][k] = solution[k];
  
        if( (number_of_generations > 0) && (q < number_of_AMS_solutions) )
        {
          out_of_range  = 1;
          shrink_factor = 2;
          while( (out_of_range == 1) && (shrink_factor > 1e-10) )
          {
            shrink_factor *= 0.5;
            out_of_range   = 0;
            for( k = 0; k < number_of_parameters; k++ )
            {
              solution_AMS[k] = solution[k] + shrink_factor*delta_AMS*distribution_multipliers[i]*(mean_vectors[i][k]-mean_vectors_previous[i][k]);
              if( !isParameterInRangeBounds( solution_AMS[k], k ) )
              {
                out_of_range = 1;
                break;
              }
            }
          }
          if( !out_of_range )
          {
            for( k = 0; k < number_of_parameters; k++ )
              populations[i][j][k] = solution_AMS[k];
          }
        }
  
//	#pragma omp task firstprivate(problem_index,i,j) shared(populations,objective_values,constraint_values) 
	{
	int info[5];
	info[0] = number_of_starts;		// start
	info[1] = 1;				// !initialize
	info[2] = number_of_generations;	// generation (step)
	info[3] = i;				// population
	info[4] = j;				// element
//        installedProblemEvaluation( problem_index, populations[i][j], &(objective_values[i][j]), &(constraint_values[i][j]) );
	torc_create(-1, installedProblemEvaluation, 5,
		1, MPI_INT, CALL_BY_COP,
		number_of_parameters, MPI_DOUBLE, CALL_BY_COP,
		1, MPI_DOUBLE, CALL_BY_RES,
		1, MPI_DOUBLE, CALL_BY_RES,
		5, MPI_INT, CALL_BY_COP,
		 &problem_index, populations[i][j], &(objective_values[i][j]), &(constraint_values[i][j]), info );
	number_of_evaluations++;

	}
  
        q++;
  
        free( solution );
      }
    }
  }
//  #pragma omp taskwait
 torc_waitall();

  }
  }

  free( solution_AMS );
}

/**
 * Computes the precision matrices required for sampling
 * the multivariate normal distribution.
 */
void computeParametersForSampling( int population_index )
{
  int i;

  if( cholesky_factors_lower_triangle[population_index] )
  {
    for( i = 0; i < number_of_parameters; i++ )
      free( cholesky_factors_lower_triangle[population_index][i] );
    free( cholesky_factors_lower_triangle[population_index] );
  }

  cholesky_factors_lower_triangle[population_index] = choleskyDecomposition( covariance_matrices[population_index], number_of_parameters );
}

/**
 * Generates and returns a single new solution by drawing
 * a single sample from a specified model.
 */
double *generateNewSolution( int population_index )
{
  short   ready;
  int     i, times_not_in_bounds;
  double *result, *z;

  times_not_in_bounds = -1;
  out_of_bounds_draws[population_index]--;

  ready = 0;
  do
  {
    times_not_in_bounds++;
    samples_drawn_from_normal[population_index]++;
    out_of_bounds_draws[population_index]++;
    if( times_not_in_bounds >= 100 )
    {
      result = (double *) Malloc( number_of_parameters*sizeof( double ) );
      for( i = 0; i < number_of_parameters; i++ )
        result[i] = lower_init_ranges[i] + (upper_init_ranges[i] - lower_init_ranges[i])*randomRealUniform01();
    }
    else
    {
      z = (double *) Malloc( number_of_parameters*sizeof( double ) );

      for( i = 0; i < number_of_parameters; i++ )
         z[i] = random1DNormalUnit();

      result = matrixVectorMultiplication( cholesky_factors_lower_triangle[population_index], z, number_of_parameters, number_of_parameters );

      for( i = 0; i < number_of_parameters; i++ )
        result[i] += mean_vectors[population_index][i];

      free( z );
    }
    
    ready = 1;
    for( i = 0; i < number_of_parameters; i++ )
    {
      if( !isParameterInRangeBounds( result[i], i ) )
      {
        ready = 0;
        break;
      }
    }
    if( !ready )
      free( result );
  }
  while( !ready );

  return( result );
}

/**
 * Adapts the distribution multipliers according to
 * the SDR-AVS mechanism.
 */
void adaptDistributionMultipliers( void )
{
  short  improvement;
  int    i;
  double st_dev_ratio;

  for( i = 0; i < number_of_populations; i++ )
  {
    if( !populations_terminated[i] )
    {
      if( (((double) out_of_bounds_draws[i])/((double) samples_drawn_from_normal[i])) > 0.9 )
        distribution_multipliers[i] *= 0.5;
  
      improvement = generationalImprovementForOnePopulation( i, &st_dev_ratio );
  
      if( improvement )
      {
        no_improvement_stretch[i] = 0;

        if( distribution_multipliers[i] < 1.0 )
          distribution_multipliers[i] = 1.0;
  
        if( st_dev_ratio > st_dev_ratio_threshold )
          distribution_multipliers[i] *= distribution_multiplier_increase;
      }
      else
      {
        if( distribution_multipliers[i] <= 1.0 )
          (no_improvement_stretch[i])++;
  
        if( (distribution_multipliers[i] > 1.0) || (no_improvement_stretch[i] >= maximum_no_improvement_stretch) )
          distribution_multipliers[i] *= distribution_multiplier_decrease;
  
        if( (no_improvement_stretch[i] < maximum_no_improvement_stretch) && (distribution_multipliers[i] < 1.0) )
          distribution_multipliers[i] = 1.0;
      }
    }
  }
}

/**
 * Determines whether an improvement is found for a specified
 * population. Returns 1 in case of an improvement, 0 otherwise.
 * The standard-deviation ratio required by the SDR-AVS
 * mechanism is computed and returned in the pointer variable.
 */
short generationalImprovementForOnePopulation( int population_index, double *st_dev_ratio )
{
  int     i, j, index_best_selected, index_best_population,
          number_of_improvements;
  double *average_parameters_of_improvements;

  /* Determine best selected solutions */
  index_best_selected = 0;
  for( i = 0; i < selection_size; i++ )
  {
    if( betterFitness( objective_values_selections[population_index][i], constraint_values_selections[population_index][i],
                       objective_values_selections[population_index][index_best_selected], constraint_values_selections[population_index][index_best_selected] ) )
      index_best_selected = i;
  }

  /* Determine best in the population and the average improvement parameters */
  average_parameters_of_improvements = (double *) Malloc( number_of_parameters*sizeof( double ) );
  for( i = 0; i < number_of_parameters; i++ )
    average_parameters_of_improvements[i] = 0.0;

  index_best_population   = 0;
  number_of_improvements  = 0;
  for( i = 0; i < population_size; i++ )
  {
    if( betterFitness( objective_values[population_index][i], constraint_values[population_index][i],
                       objective_values[population_index][index_best_population], constraint_values[population_index][index_best_population] ) )
      index_best_population = i;

    if( betterFitness( objective_values[population_index][i], constraint_values[population_index][i],
                       objective_values_selections[population_index][index_best_selected], constraint_values_selections[population_index][index_best_selected] ) )
    {
      number_of_improvements++;
      for( j = 0; j < number_of_parameters; j++ )
        average_parameters_of_improvements[j] += populations[population_index][i][j];
    }
  }

  /* Determine st.dev. ratio */
  *st_dev_ratio = 0.0;
  if( number_of_improvements > 0 )
  {
    for( i = 0; i < number_of_parameters; i++ )
      average_parameters_of_improvements[i] /= (double) number_of_improvements;

    *st_dev_ratio = getStDevRatio( population_index, average_parameters_of_improvements );
  }

  free( average_parameters_of_improvements );

  if( fabs( objective_values_selections[population_index][index_best_selected] - objective_values[population_index][index_best_population] ) == 0.0 )
    return( 0 );

  return( 1 );
}

/**
 * Returns 1 if x is better than y, 0 otherwise.
 * x is not better than y unless:
 * - x and y are both infeasible and x has a smaller sum of constraint violations, or
 * - x is feasible and y is not, or
 * - x and y are both feasible and x has a smaller objective value than y
 */
short betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
  short result;
  
  result = 0;

  if( constraint_value_x > 0 ) /* x is infeasible */
  {
    if( constraint_value_y > 0 ) /* Both are infeasible */
    {
      if( constraint_value_x < constraint_value_y )
       result = 1;
    }
  }
  else /* x is feasible */
  {
    if( constraint_value_y > 0 ) /* x is feasible and y is not */
      result = 1;
    else /* Both are feasible */
    {
      if( objective_value_x < objective_value_y )
        result = 1;
    }
  }

  return( result );
}

/**
 * Computes and returns the standard-deviation-ratio
 * of a given point for a given model.
 */
double getStDevRatio( int population_index, double *parameters )
{
  int      i;
  double **inverse, result, *x_min_mu, *z;

  inverse = matrixLowerTriangularInverse( cholesky_factors_lower_triangle[population_index], number_of_parameters );

  x_min_mu = (double *) Malloc( number_of_parameters*sizeof( double ) );
  
  for( i = 0; i < number_of_parameters; i++ )
    x_min_mu[i] = parameters[i]-mean_vectors[population_index][i];

  z = matrixVectorMultiplication( inverse, x_min_mu, number_of_parameters, number_of_parameters );

  result = 0.0;
  for( i = 0; i < number_of_parameters; i++ )
  {
    if( fabs( z[i] ) > result )
      result = fabs( z[i] );
  }

  free( z );
  free( x_min_mu );
  for( i = 0; i < number_of_parameters; i++ )
    free( inverse[i] );
  free( inverse );

  return( result );
}

/**
 * Determines the best solution of the current run
 * and compares it to the best solution found in all runs.
 * Then stores the very best solution found so far.
 */
void determineBestSolutionSoFar( void )
{
  int i, population_of_best, index_of_best;

  determineBestSolutionInCurrentPopulations( &population_of_best, &index_of_best );

  if( number_of_starts == 1 ||
      betterFitness( objective_values[population_of_best][index_of_best],
                     constraint_values[population_of_best][index_of_best],
                     best_so_far_objective_value,
                     best_so_far_constraint_value ) )
  {
    best_so_far_objective_value  = objective_values[population_of_best][index_of_best];
    best_so_far_constraint_value = constraint_values[population_of_best][index_of_best];
    for( i = 0; i < number_of_parameters; i++ )
      best_so_far_solution[i] = populations[population_of_best][index_of_best][i];
  }
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Ezilaitini -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitini( void )
{
  ezilaitiniMemory();

  ezilaitiniDistributionMultipliers();
  
  ezilaitiniObjectiveRotationMatrix();
}

/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitiniMemory( void )
{
  int i, j;

  for( i = 0; i < number_of_populations; i++ )
  {
    for( j = 0; j < population_size; j++ )
      free( populations[i][j] );
    free( populations[i] );
    
    free( objective_values[i] );
    
    free( constraint_values[i] );

    free( ranks[i] );

    for( j = 0; j < selection_size; j++ )
      free( selections[i][j] );
    free( selections[i] );

    free( objective_values_selections[i] );

    free( constraint_values_selections[i] );

    free( mean_vectors[i] );

    free( mean_vectors_previous[i] );

    for( j = 0; j < number_of_parameters; j++ )
      free( covariance_matrices[i][j] );
    free( covariance_matrices[i] );

    if( cholesky_factors_lower_triangle[i] )
    {
      for( j = 0; j < number_of_parameters; j++ )
        free( cholesky_factors_lower_triangle[i][j] );
      free( cholesky_factors_lower_triangle[i] );
    }
  }

  free( covariance_matrices );
  free( cholesky_factors_lower_triangle );
  free( lower_range_bounds );
  free( upper_range_bounds );
  free( lower_init_ranges );
  free( upper_init_ranges );
  free( populations_terminated );
  free( no_improvement_stretch );
  free( populations );
  free( objective_values );
  free( constraint_values );
  free( ranks );
  free( selections );
  free( objective_values_selections );
  free( constraint_values_selections );
  free( mean_vectors );
  free( mean_vectors_previous );
}

/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitiniDistributionMultipliers( void )
{
  free( distribution_multipliers );
  free( samples_drawn_from_normal );
  free( out_of_bounds_draws );
}

/**
 * Undoes initialization procedure by freeing up memory.
 */
void ezilaitiniObjectiveRotationMatrix( void )
{
  int i;

  if( rotation_angle == 0.0 )
    return;

  for( i = 0; i < number_of_parameters; i++ )
    free( rotation_matrix[i] );
  free( rotation_matrix );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Runs the IDEA with exponentially
 * increasing population size and
 * number of parallel populations.
 */
void run( void )
{
  double t0, t1;
  int population_size_base;

  number_of_starts = 0;

  best_so_far_solution = (double *) Malloc( number_of_parameters*sizeof( double ) );

  distribution_multiplier_decrease = 0.9;
  st_dev_ratio_threshold           = 1.0;
  tau                              = 0.35;
  maximum_no_improvement_stretch   = 25+number_of_parameters;
  population_size_base             = (int) (17.0 + 3.0*pow((double) number_of_parameters,1.5));

  do
  {
    if( number_of_starts % 2 == 0 )
    {
      population_size       = (1+number_of_starts/2)*population_size_base;
      number_of_populations = (int) pow(2.0,number_of_starts/2);

      if( number_of_populations > maximum_number_of_populations )
      {
        population_size       = (population_size*number_of_populations)/maximum_number_of_populations;
        number_of_populations = maximum_number_of_populations;
      }
    }
    else
    {
      population_size       = ((int) pow(2.0,1+(number_of_starts/2)))*population_size_base;
      number_of_populations = 1;
    }

    t0 = torc_gettime();
    runOnce();
    t1 = torc_gettime();
    printf("step %d in %.3f seconds\n", number_of_starts, t1-t0); 
  }
  while( !checkTerminationCondition() );

  free( best_so_far_solution );
}

/**
 * Runs the IDEA once with fixed population size
 * and number of parallel populations.
 */
void runOnce( void )
{
  number_of_starts++;

  initialize();

  if( print_verbose_overview && (number_of_starts == 1) )
    printVerboseOverview();

  while( !checkTerminationConditionForRunOnce() )
  {
    if( write_generational_statistics )
      writeGenerationalStatistics();

    if( write_generational_solutions )
      writeGenerationalSolutions( 0 );

    makeSelections();
    
    makePopulations();

    number_of_generations++;
  }

  determineBestSolutionSoFar();

  if( checkTerminationCondition() )
  {
    writeGenerationalStatistics();

    writeGenerationalSolutions( 1 );
  }

  ezilaitini();
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Main -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * The main function:
 * - interpret parameters on the command line
 * - run the algorithm with the interpreted parameters
 */
int main( int argc, char **argv )
{
  interpretCommandLine( argc, argv );

  torc_register_task(installedProblemEvaluation);
  torc_init(argc, argv, MODE_MW);

  double t0 = torc_gettime();
  run();
  double t1 = torc_gettime();
  printf("Total elapsed time = %.3f seconds\n", t1-t0);

  torc_finalize();
  return( 0 );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
