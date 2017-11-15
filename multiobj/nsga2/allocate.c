/* Memory allocation and deallocation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to allocate memory to a population */
void allocate_memory_pop (population *pop, int size)
{
    int i;
    pop->ind = (individual *)calloc(1, size*sizeof(individual));
    for (i=0; i<size; i++)
    {
        allocate_memory_ind (&(pop->ind[i]));
    }
    return;
}

/* Function to allocate memory to an individual */
void allocate_memory_ind (individual *ind)
{
    int j;
    if (nreal != 0)
    {
        ind->xreal = (double *)calloc(1, nreal*sizeof(double));
    }
    if (nbin != 0)
    {
        ind->xbin = (double *)calloc(1, nbin*sizeof(double));
        ind->gene = (int **)calloc(1, nbin*sizeof(int *));
        for (j=0; j<nbin; j++)
        {
            ind->gene[j] = (int *)calloc(1, nbits[j]*sizeof(int));
        }
    }
    ind->obj = (double *)calloc(1, nobj*sizeof(double));
    if (ncon != 0)
    {
        ind->constr = (double *)calloc(1, ncon*sizeof(double));
    }
    return;
}

/* Function to deallocate memory to a population */
void deallocate_memory_pop (population *pop, int size)
{
    int i;
    for (i=0; i<size; i++)
    {
        deallocate_memory_ind (&(pop->ind[i]));
    }
    free (pop->ind);
    return;
}

/* Function to deallocate memory to an individual */
void deallocate_memory_ind (individual *ind)
{
    int j;
    if (nreal != 0)
    {
        free(ind->xreal);
    }
    if (nbin != 0)
    {
        for (j=0; j<nbin; j++)
        {
            free(ind->gene[j]);
        }
        free(ind->xbin);
        free(ind->gene);
    }
    free(ind->obj);
    if (ncon != 0)
    {
        free(ind->constr);
    }
    return;
}
