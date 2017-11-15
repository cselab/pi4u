/* Data initializtion routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to initialize a population randomly */
void initialize_pop (population *pop)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        initialize_ind (&(pop->ind[i]));
    }
    return;
}

/* Function to initialize an individual randomly */
void initialize_ind (individual *ind)
{
    int j, k;
    if (nreal!=0)
    {
        for (j=0; j<nreal; j++)
        {
            ind->xreal[j] = rndreal (min_realvar[j], max_realvar[j]);
        }
    }
    if (nbin!=0)
    {
        for (j=0; j<nbin; j++)
        {
            for (k=0; k<nbits[j]; k++)
            {
                if (randomperc() <= 0.5)
                {
                    ind->gene[j][k] = 0;
                }
                else
                {
                    ind->gene[j][k] = 1;
                }
            }
        }
    }
    return;
}

void initialize_pop_fp (population *pop, FILE *fp)
{
    int i;
    for (i=0; i<popsize; i++)
    {
        initialize_ind_fp (&(pop->ind[i]), fp, i);
    }
    return;
}

/* Function to initialize an individual randomly */
void initialize_ind_fp (individual *ind, FILE *fp, int i)
{
    int j, k;
    if (nreal!=0)
    {
        for (j=0; j<nreal; j++)
        {
#if 0
            ind->xreal[j] = rndreal (min_realvar[j], max_realvar[j]);
#else
            fscanf(fp, "%lf", &ind->xreal[j]);
	    printf("read %d %d : %lf\n", i, j, ind->xreal[j]);
#endif
        }
    }
    if (nbin!=0)
    {
        for (j=0; j<nbin; j++)
        {
            for (k=0; k<nbits[j]; k++)
            {
                if (randomperc() <= 0.5)
                {
                    ind->gene[j][k] = 0;
                }
                else
                {
                    ind->gene[j][k] = 1;
                }
            }
        }
    }
    return;
}
