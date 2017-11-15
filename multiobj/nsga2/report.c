/* Routines for storing population data into files */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"

/* Function to print the information of a population in a file */
void report_pop (population *pop, FILE *fpt)
{
    int i, j, k;
    for (i=0; i<popsize; i++)
    {
        for (j=0; j<nobj; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
        }
        if (ncon!=0)
        {
            for (j=0; j<ncon; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
            }
        }
        if (nreal!=0)
        {
            for (j=0; j<nreal; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
            }
        }
        if (nbin!=0)
        {
            for (j=0; j<nbin; j++)
            {
                for (k=0; k<nbits[j]; k++)
                {
                    fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                }
            }
        }
        fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
        fprintf(fpt,"%d\t",pop->ind[i].rank);
        fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
    }
    return;
}

void dump_pop_txt (population *pop, FILE *fpt)
{
    for (int i=0; i<popsize; i++)
    {
        if (nreal!=0)
        {
            for (int j=0; j<nreal; j++)
            {
                fprintf(fpt,"%.10f ",pop->ind[i].xreal[j]);
            }
        }

        for (int j=0; j<nobj; j++)
        {
            fprintf(fpt,"%.10f ",pop->ind[i].obj[j]);
        }
#if 0
        if (ncon!=0)
        {
            for (int j=0; j<ncon; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
            }
        }
#endif
#if 0
        if (nbin!=0)
        {
            for (int j=0; j<nbin; j++)
            {
                for (int k=0; k<nbits[j]; k++)
                {
                    fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                }
            }
        }
#endif
        fprintf(fpt,"%.10f ",pop->ind[i].constr_violation);
#if 0
        fprintf(fpt,"%d ",pop->ind[i].rank);
        fprintf(fpt,"%f ",pop->ind[i].crowd_dist);
#endif
	fprintf(fpt, "\n");
    }
    return;
}

void dump_pop (population *pop, FILE *fpt)
{
    for (int i=0; i<popsize; i++)
    {
        for (int j=0; j<nobj; j++)
        {
            fwrite(&pop->ind[i].obj[j], sizeof(double), 1, fpt);
        }
        if (nreal!=0)
        {
            for (int j=0; j<nreal; j++)
            {
		fwrite(&pop->ind[i].xreal[j], sizeof(double), 1, fpt);
            }
        }
	fwrite(&pop->ind[i].constr_violation, sizeof(double), 1, fpt);
    }
    return;
}

void load_pop (population *pop, FILE *fpt)
{
    for (int i=0; i<popsize; i++)
    {
        for (int j=0; j<nobj; j++)
        {
            fread(&pop->ind[i].obj[j], sizeof(double), 1, fpt);
        }
        if (nreal!=0)
        {
            for (int j=0; j<nreal; j++)
            {
		fread(&pop->ind[i].xreal[j], sizeof(double), 1, fpt);
            }
        }
	fread(&pop->ind[i].constr_violation, sizeof(double), 1, fpt);
    }
    return;
}

void load_pop_txt (population *pop, FILE *fpt)
{
    for (int i=0; i<popsize; i++)
    {
        if (nreal!=0)
        {
            for (int j=0; j<nreal; j++)
            {
                fscanf(fpt,"%lf",&pop->ind[i].xreal[j]);
            }
        }
        for (int j=0; j<nobj; j++)
        {
            fscanf(fpt,"%lf",&pop->ind[i].obj[j]);
        }
#if 0
        if (ncon!=0)
        {
            for (int j=0; j<ncon; j++)
            {
                fscanf(fpt,"%lf",&pop->ind[i].constr[j]);
            }
        }
#endif
#if 0
        if (nbin!=0)
        {
            for (int j=0; j<nbin; j++)
            {
                for (int k=0; k<nbits[j]; k++)
                {
                    fscanf(fpt,"%d",&pop->ind[i].gene[j][k]);
                }
            }
        }
#endif
        fscanf(fpt,"%lf",&pop->ind[i].constr_violation);
#if 0
        fscanf(fpt,"%d",&pop->ind[i].rank);
        fscanf(fpt,"%lf",&pop->ind[i].crowd_dist);
#endif
    }
    return;
}


/* Function to read the information of a population from a file */
void load_pop0 (population *pop, FILE *fpt)
{
    for (int i=0; i<popsize; i++)
    {
        for (int j=0; j<nobj; j++)
        {
            fscanf(fpt,"%le",&pop->ind[i].obj[j]);
        }
        if (ncon!=0)
        {
            for (int j=0; j<ncon; j++)
            {
                fscanf(fpt,"%le",&pop->ind[i].constr[j]);
            }
        }
        if (nreal!=0)
        {
            for (int j=0; j<nreal; j++)
            {
                fscanf(fpt,"%le",&pop->ind[i].xreal[j]);
            }
        }
        if (nbin!=0)
        {
            for (int j=0; j<nbin; j++)
            {
                for (int k=0; k<nbits[j]; k++)
                {
                    fscanf(fpt,"%d",&pop->ind[i].gene[j][k]);
                }
            }
        }
        fscanf(fpt,"%le",&pop->ind[i].constr_violation);
        fscanf(fpt,"%d",&pop->ind[i].rank);
        fscanf(fpt,"%le",&pop->ind[i].crowd_dist);
    }
    return;
}

/* Function to print the information of feasible and non-dominated population in a file */
void report_feasible (population *pop, FILE *fpt)
{
    for (int i=0; i<popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            for (int j=0; j<nobj; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
            }
            if (ncon!=0)
            {
                for (int j=0; j<ncon; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
                }
            }
            if (nreal!=0)
            {
                for (int j=0; j<nreal; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
                }
            }
            if (nbin!=0)
            {
                for (int j=0; j<nbin; j++)
                {
                    for (int k=0; k<nbits[j]; k++)
                    {
                        fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                    }
                }
            }
            fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\t",pop->ind[i].rank);
            fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
        }
    }
    return;
}
