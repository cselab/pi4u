/* NSGA-II routine (implementation of the 'main' function) */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <unistd.h>
# include <string.h>
# include "global.h"
# include "rand.h"

#if defined(_USE_TORC_)
# include <torc.h>
#endif

#if !defined(_USE_TORC_)
#include <sys/time.h>
static double torc_gettime()
{
        struct timeval t;
        gettimeofday(&t, NULL);
        return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}
#endif


/* JOBMAXTIME is used for clean restarts, at the moment works only if 
 * first two generations are executed in same job */
#define JOBMAXTIME 24*60*60	// 24 hours in seconds
				
int nreal;
int nbin;
int nobj;
int ncon;
int popsize;
double pcross_real;
double pcross_bin;
double pmut_real;
double pmut_bin;
double eta_c;
double eta_m;
int ngen;
int nbinmut;
int nrealmut;
int nbincross;
int nrealcross;
int *nbits;
double *min_realvar;
double *max_realvar;
double *min_binvar;
double *max_binvar;
int bitlength;
int choice;
int obj1;
int obj2;
int obj3;
int angle1;
int angle2;

int main (int argc, char **argv)
{
    int checkpoint_restart = 0;

    double gt0, gt1, gt2;	// used for clean restart mechanism
    int i;
    FILE *fpt1;
    FILE *fpt2;
    FILE *fpt3;
    FILE *fpt4;
    FILE *fpt5;
    FILE *gp = NULL;
    population *parent_pop;
    population *child_pop;
    population *mixed_pop;

#if defined(_USE_TORC_)
	torc_register_task(test_problem_v2);
	torc_init(argc, argv, MODE_MW);
#endif

	gt0 = torc_gettime();	// time at start of execution

    if (argc<2)
    {
        printf("\n Usage ./nsga2r <config.in> \n");
        exit(1);
    }

    FILE *fpin = fopen(argv[1], "r");
    if (fpin == NULL)
    {
        printf("Failed to open configuration file %s. Exiting!\n", argv[1]);
        return 1;
    }

    printf("\n Enter the random seed (0-1): ");
    fscanf(fpin, "%lf", &seed);
    if (fpin != stdout) printf("%lf\n", seed);
    if (seed<=0.0 || seed>=1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }
    fpt1 = fopen("initial_pop.out","w");
    fpt2 = fopen("final_pop.out","w");
    fpt3 = fopen("best_pop.out","w");
    fpt4 = fopen("all_pop.out","w");
    fpt5 = fopen("params.out","w");
    fprintf(fpt1,"# This file contains the data of initial population\n");
    fprintf(fpt2,"# This file contains the data of final population\n");
    fprintf(fpt3,"# This file contains the data of final feasible population (if found)\n");
    fprintf(fpt4,"# This file contains the data of all generations\n");
    fprintf(fpt5,"# This file contains information about inputs as read by the program\n");
    printf("\n Enter the problem relevant and algorithm relevant parameters ... ");

    printf("\n Enter the population size (a multiple of 4) : ");
    fscanf(fpin, "%d", &popsize);
    if (fpin != stdout) printf("%d\n", popsize);
    if (popsize<4 || (popsize%4)!= 0)
    {
        printf("\n population size read is : %d",popsize);
        printf("\n Wrong population size entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of generations : ");
    fscanf(fpin, "%d", &ngen);
    printf("%d\n", ngen);
    if (ngen<1)
    {
        printf("\n number of generations read is : %d",ngen);
        printf("\n Wrong nuber of generations entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of objectives : ");
    fscanf(fpin, "%d", &nobj);
    if (fpin != stdout) printf("%d\n", nobj);
    if (nobj<1)
    {
        printf("\n number of objectives entered is : %d",nobj);
        printf("\n Wrong number of objectives entered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of constraints : ");
    fscanf(fpin, "%d", &ncon);
    if (fpin != stdout) printf("%d\n", ncon);
    if (ncon<0)
    {
        printf("\n number of constraints entered is : %d",ncon);
        printf("\n Wrong number of constraints enetered, hence exiting \n");
        exit (1);
    }
    printf("\n Enter the number of real variables : ");
    fscanf(fpin, "%d", &nreal);
    if (fpin != stdout) printf("%d\n", nreal);
    if (nreal<0)
    {
        printf("\n number of real variables entered is : %d",nreal);
        printf("\n Wrong number of variables entered, hence exiting \n");
        exit (1);
    }
    if (nreal != 0)
    {
        min_realvar = (double *)malloc(nreal*sizeof(double));
        max_realvar = (double *)malloc(nreal*sizeof(double));
        for (i=0; i<nreal; i++)
        {
            printf ("\n Enter the lower limit of real variable %d : ",i+1);
            fscanf(fpin, "%lf", &min_realvar[i]);
            if (fpin != stdout) printf(" %lf\n", min_realvar[i]);
            printf ("\n Enter the upper limit of real variable %d : ",i+1);
            fscanf(fpin, "%lf",&max_realvar[i]);
            if (fpin != stdout) printf(" %lf\n", max_realvar[i]);
            if (max_realvar[i] <= min_realvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of real variable (0.6-1.0) : ");
        fscanf(fpin, "%lf",&pcross_real);
        if (fpin != stdout) printf("%lf\n", pcross_real);
        if (pcross_real<0.0 || pcross_real>1.0)
        {
            printf("\n Probability of crossover entered is : %e",pcross_real);
            printf("\n Entered value of probability of crossover of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probablity of mutation of real variables (1/nreal) : ");
        fscanf(fpin, "%lf",&pmut_real);
        if (fpin != stdout) printf("%lf\n", pmut_real);
        if (pmut_real<0.0 || pmut_real>1.0)
        {
            printf("\n Probability of mutation entered is : %e",pmut_real);
            printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for crossover (5-20): ");
        fscanf(fpin, "%lf", &eta_c);
        if (fpin != stdout) printf("%lf\n", eta_c);
        if (eta_c<=0)
        {
            printf("\n The value entered is : %e",eta_c);
            printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for mutation (5-50): ");
        fscanf(fpin, "%lf", &eta_m);
        if (fpin != stdout) printf("%lf\n", eta_m);
        if (eta_m<=0)
        {
            printf("\n The value entered is : %e",eta_m);
            printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
            exit (1);
        }
    }
    printf("\n Enter the number of binary variables : ");
    fscanf(fpin, "%d", &nbin);
    if (fpin != stdout) printf("%d\n", nbin);
    if (nbin<0)
    {
        printf ("\n number of binary variables entered is : %d",nbin);
        printf ("\n Wrong number of binary variables entered, hence exiting \n");
        exit(1);
    }
    if (nbin != 0)
    {
#if 1
        printf ("\n number of binary variables entered is : %d",nbin);
        printf ("\n This option has been currently disabled in the parallel version\n");
        exit(1);
	
#endif
        nbits = (int *)malloc(nbin*sizeof(int));
        min_binvar = (double *)malloc(nbin*sizeof(double));
        max_binvar = (double *)malloc(nbin*sizeof(double));
        for (i=0; i<nbin; i++)
        {
            printf ("\n Enter the number of bits for binary variable %d : ",i+1);
            fscanf(fpin, "%d", &nbits[i]);
            if (fpin != stdout) printf("%d\n", nbits[i]);
            if (nbits[i] < 1)
            {
                printf("\n Wrong number of bits for binary variable entered, hence exiting");
                exit(1);
            }
            printf ("\n Enter the lower limit of binary variable %d : ",i+1);
            fscanf(fpin, "%lf", &min_binvar[i]);
            if (fpin != stdout) printf("%lf\n", min_binvar[i]);
            printf ("\n Enter the upper limit of binary variable %d : ",i+1);
            fscanf(fpin, "%lf", &max_binvar[i]);
            if (fpin != stdout) printf("%lf\n", max_binvar[i]);
            if (max_binvar[i] <= min_binvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
        fscanf(fpin, "%lf", &pcross_bin);
        if (fpin != stdout) printf("%lf\n", pcross_bin);
        if (pcross_bin<0.0 || pcross_bin>1.0)
        {
            printf("\n Probability of crossover entered is : %e",pcross_bin);
            printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probability of mutation of binary variables (1/nbits): ");
        fscanf(fpin, "%lf", &pmut_bin);
        if (fpin != stdout) printf("%lf\n", pmut_bin);
        if (pmut_bin<0.0 || pmut_bin>1.0)
        {
            printf("\n Probability of mutation entered is : %e",pmut_bin);
            printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
    }
    if (nreal==0 && nbin==0)
    {
        printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
        exit(1);
    }
    choice=0;
    printf("\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) : ");
    fscanf(fpin, "%d",&choice);
    if (fpin != stdout) printf("%d\n", choice);
    if (choice!=0 && choice!=1)
    {
        printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
        exit(1);
    }
    if (choice==1)
    {
        gp = popen(GNUPLOT_COMMAND,"w");
        if (gp==NULL)
        {
            printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
            printf("\n Edit the string to suit your system configuration and rerun the program\n");
            exit(1);
        }
        if (nobj==2)
        {
            printf("\n Enter the objective for X axis display : ");
            fscanf(fpin, "%d", &obj1);
            if (fpin != stdout) printf("%d\n", obj1); 
            if (obj1<1 || obj1>nobj)
            {
                printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                exit(1);
            }
            printf("\n Enter the objective for Y axis display : ");
            fscanf(fpin, "%d",&obj2);
            if (fpin != stdout) printf("%d\n", obj2); 
            if (obj2<1 || obj2>nobj)
            {
                printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                exit(1);
            }
            obj3 = -1;
        }
        else
        {
            printf("\n #obj > 2, 2D display or a 3D display ?, enter 2 for 2D and 3 for 3D :");
            fscanf(fpin, "%d",&choice);
            if (choice!=2 && choice!=3)
            {
                printf("\n Entered the wrong choice, hence exiting, choice entered was %d\n",choice);
                exit(1);
            }
            if (choice==2)
            {
                printf("\n Enter the objective for X axis display : ");
                fscanf(fpin, "%d",&obj1);
                if (fpin != stdout) printf("%d\n", obj1); 
                if (obj1<1 || obj1>nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                fscanf(fpin, "%d",&obj2);
                if (fpin != stdout) printf("%d\n", obj2); 
                if (obj2<1 || obj2>nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                    exit(1);
                }
                obj3 = -1;
            }
            else
            {
                printf("\n Enter the objective for X axis display : ");
                fscanf(fpin, "%d",&obj1);
                if (fpin != stdout) printf("%d\n", obj1); 
                if (obj1<1 || obj1>nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                fscanf(fpin, "%d",&obj2);
                if (fpin != stdout) printf("%d\n", obj2); 
                if (obj2<1 || obj2>nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",obj2);
                    exit(1);
                }
                printf("\n Enter the objective for Z axis display : ");
                fscanf(fpin, "%d",&obj3);
                if (fpin != stdout) printf("%d\n", obj3);
                if (obj3<1 || obj3>nobj)
                {
                    printf("\n Wrong value of Z objective entered, value entered was %d\n",obj3);
                    exit(1);
                }
                printf("\n You have chosen 3D display, hence location of eye required \n");
                printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
                fscanf(fpin, "%d",&angle1);
                if (fpin != stdout) printf("%d\n", angle1);
                if (angle1<0 || angle1>180)
                {
                    printf("\n Wrong value for first angle entered, hence exiting \n");
                    exit(1);
                }
                printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
                fscanf(fpin, "%d",&angle2);
                if (fpin != stdout) printf("%d\n", angle2);
                if (angle2<0 || angle2>360)
                {
                    printf("\n Wrong value for second angle entered, hence exiting \n");
                    exit(1);
                }
            }
        }
    }
    printf("\n Input data successfully entered, now performing initialization \n");
    fprintf(fpt5,"\n Population size = %d",popsize);
    fprintf(fpt5,"\n Number of generations = %d",ngen);
    fprintf(fpt5,"\n Number of objective functions = %d",nobj);
    fprintf(fpt5,"\n Number of constraints = %d",ncon);
    fprintf(fpt5,"\n Number of real variables = %d",nreal);
    if (nreal!=0)
    {
        for (i=0; i<nreal; i++)
        {
            fprintf(fpt5,"\n Lower limit of real variable %d = %e",i+1,min_realvar[i]);
            fprintf(fpt5,"\n Upper limit of real variable %d = %e",i+1,max_realvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of real variable = %e",pcross_real);
        fprintf(fpt5,"\n Probability of mutation of real variable = %e",pmut_real);
        fprintf(fpt5,"\n Distribution index for crossover = %e",eta_c);
        fprintf(fpt5,"\n Distribution index for mutation = %e",eta_m);
    }
    fprintf(fpt5,"\n Number of binary variables = %d",nbin);
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            fprintf(fpt5,"\n Number of bits for binary variable %d = %d",i+1,nbits[i]);
            fprintf(fpt5,"\n Lower limit of binary variable %d = %e",i+1,min_binvar[i]);
            fprintf(fpt5,"\n Upper limit of binary variable %d = %e",i+1,max_binvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of binary variable = %e",pcross_bin);
        fprintf(fpt5,"\n Probability of mutation of binary variable = %e",pmut_bin);
    }
    fprintf(fpt5,"\n Seed for random number generator = %e",seed);
    bitlength = 0;
    if (nbin!=0)
    {
        for (i=0; i<nbin; i++)
        {
            bitlength += nbits[i];
        }
    }
    fprintf(fpt1,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt2,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt3,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    fprintf(fpt4,"# of objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nobj,ncon,nreal,bitlength);
    nbinmut = 0;
    nrealmut = 0;
    nbincross = 0;
    nrealcross = 0;
    parent_pop = (population *)calloc(1, sizeof(population));
    child_pop = (population *)calloc(1, sizeof(population));
    mixed_pop = (population *)calloc(1, sizeof(population));
    allocate_memory_pop (parent_pop, popsize);
    allocate_memory_pop (child_pop, popsize);
    allocate_memory_pop (mixed_pop, 2*popsize);
    randomize();

	char user_initial_pop[256];
	strcpy(user_initial_pop, "user_initial_pop.txt");
	FILE *fp = fopen(user_initial_pop, "r");
	if (fp == NULL)
	{
		initialize_pop (parent_pop);
        
	}
	else
	{
		printf("Initializing first population from user file...\n");
		if (checkpoint_restart == 1)
		{
			printf("Setting checkpoint_restart flag to 0!\n");
			checkpoint_restart = 0;
		}
		initialize_pop_fp (parent_pop, fp);
		fclose(fp);
	}

    	printf("\n Initialization done, now performing first generation\n");
    	decode_pop(parent_pop);

	if (checkpoint_restart)
	{
		FILE *fp;
		char fname[256];
		sprintf(fname, "eval_db_%03d.txt", 1);
#if BINARY
		fp = fopen(fname, "rb");
#else
		fp = fopen(fname, "r");
#endif
		if (fp != NULL)
		{
			printf("loading pop %d\n", 1);
#if BINARY
			load_pop(parent_pop, fp);
#else
			load_pop_txt(parent_pop, fp);
#endif
			fclose(fp);
		}
		else
		{
			printf("done loading\n");
			checkpoint_restart = 0;
		}
	}

	if (!checkpoint_restart)
	{
		printf("evaluating pop %d\n", 1);
		
		gt1 = torc_gettime();
		evaluate_pop (parent_pop, 1);
		gt2 = torc_gettime();
		/*could dump here*/
	}

    	assign_rank_and_crowding_distance (parent_pop);

#if IODUMP
	if (!checkpoint_restart)
	{
		FILE *fp;
		char fname[256];
		sprintf(fname, "eval_db_%03d.txt", 1);
		printf("dumping pop %d in eval_db \n", 1);
#if BINARY
		fp = fopen(fname, "wb");
		dump_pop(parent_pop, fp);
#else
		fp = fopen(fname, "w");
		dump_pop_txt(parent_pop, fp);
#endif
		fclose(fp);
	}

#endif
	// avoid repeated dumping
	if (!checkpoint_restart) 
	{
    		report_pop (parent_pop, fpt1);
    		fprintf(fpt4,"# gen = 1\n");
    		report_pop(parent_pop,fpt4);
#if IODUMP
		{
		FILE *fp;
		char fname[256];
		sprintf(fname, "curgen_db_%03d.txt", 1);
		printf("dumping pop %d in curgend_db \n", 1);
		fp = fopen(fname, "w");
		dump_pop_txt(parent_pop, fp);
		fclose(fp);
		}
	}
#endif
    printf("\n gen = 1 done \n");
    fflush(stdout);
    if (choice!=0)    onthefly_display (parent_pop,gp,1);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    sleep(1);

    for (i=2; i<=ngen; i++)
    {
	selection (parent_pop, child_pop);
	mutation_pop (child_pop);
	decode_pop(child_pop);

#if 1
	if (checkpoint_restart)
	{
		FILE *fp;
		char fname[256];
		sprintf(fname, "eval_db_%03d.txt", i);
#if BINARY
		fp = fopen(fname, "rb");
#else
		fp = fopen(fname, "r");
#endif
		if (fp != NULL)
		{
			printf("loading pop %d\n", i);
#if BINARY
			load_pop(child_pop, fp);
#else
			load_pop_txt(child_pop, fp);
#endif
			fclose(fp);
		}
		else
		{
			checkpoint_restart = 0;
		}
	}

	if (!checkpoint_restart)
	{
		printf("evaluating pop %d\n", i);
		gt1 = torc_gettime();
		evaluate_pop (child_pop, i);
		gt2 = torc_gettime();

#if IODUMP
		{
		FILE *fp;
		char fname[256];
		sprintf(fname, "eval_db_%03d.txt", i);
		printf("dumping pop %d\n", i);
#if BINARY
		fp = fopen(fname, "wb");
		dump_pop(child_pop, fp);
#else
		fp = fopen(fname, "w");
		dump_pop_txt(child_pop, fp);
#endif
		fclose(fp);
		}
#endif

	}
#else
        evaluate_pop(child_pop, i);
#endif

        merge (parent_pop, child_pop, mixed_pop);
        fill_nondominated_sort (mixed_pop, parent_pop);
 	
	// the following if statement has to stay after merge and fillnds!!
	if (!checkpoint_restart) 
	{
		/* Comment following four lines if information for all
		generations is not desired, it will speed up the execution */
		fprintf(fpt4,"# gen = %d\n",i);
		report_pop(parent_pop,fpt4);
        	fflush(fpt4);
        	if (choice!=0)    onthefly_display (parent_pop,gp,i);
#if IODUMP
		{
		FILE *fp;
		char fname[256];
		sprintf(fname, "curgen_db_%03d.txt", i);
		printf("dumping pop %d\n", i);
		fp = fopen(fname, "w");
		dump_pop_txt(parent_pop, fp);
		fclose(fp);
		}
	}
#endif
       
        printf("\n gen = %d done\n",i);

	// check whether there is enough time to perform another generation
#if (JOBMAXTIME > 0)
	{
	double lastgen_time = gt2 - gt1;
	static long maxt = -1;

	long runt, remt;
	if (maxt == -1) {
		maxt = JOBMAXTIME;
		printf("job maxtime = %ld\n", maxt);
	}

	runt = torc_gettime() - gt0;
	remt = maxt - runt;
	printf("job runtime = %ld\n", runt);
	printf("job remtime = %ld\n", remt);

	if ((lastgen_time * 1.5) > remt) {
		printf("no more available time, exiting...\n");
		// do something
		break;		
	}

	}
#endif	
    }
    printf("\n Generations finished, now reporting solutions");
    report_pop(parent_pop,fpt2);
    report_feasible(parent_pop,fpt3);
    if (nreal!=0)
    {
        fprintf(fpt5,"\n Number of crossover of real variable = %d",nrealcross);
        fprintf(fpt5,"\n Number of mutation of real variable = %d",nrealmut);
    }
    if (nbin!=0)
    {
        fprintf(fpt5,"\n Number of crossover of binary variable = %d",nbincross);
        fprintf(fpt5,"\n Number of mutation of binary variable = %d",nbinmut);
    }
    fflush(stdout);
    fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    fflush(fpt4);
    fflush(fpt5);
    fclose(fpt1);
    fclose(fpt2);
    fclose(fpt3);
    fclose(fpt4);
    fclose(fpt5);
    if (choice!=0)
    {
        pclose(gp);
    }
    if (nreal!=0)
    {
        free (min_realvar);
        free (max_realvar);
    }
    if (nbin!=0)
    {
        free (min_binvar);
        free (max_binvar);
        free (nbits);
    }
    deallocate_memory_pop (parent_pop, popsize);
    deallocate_memory_pop (child_pop, popsize);
    deallocate_memory_pop (mixed_pop, 2*popsize);
    free (parent_pop);
    free (child_pop);
    free (mixed_pop);
    printf("\n Routine successfully exited \n");


#if defined(_USE_TORC_)
	torc_finalize();
#endif
    return (0);
}
