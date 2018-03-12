#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <math.h>
#include <pndlc.h>
#include <torc.h>

static double my_gettime()
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}


double F(double *TP, int *pn)
{
    double gres = 0;
    int i, n = *pn;

    printf("FUN running for %.8lf %.8lf %.8lf\n", TP[0], TP[1], TP[2]);
    //sleep(1);
    gres = 0;
    for (i = 0; i<n; i++)    // x^2+y^2+z^2
        gres += TP[i]*TP[i]*TP[i];

    gres += TP[0]*TP[1]*TP[2];    // x*y

    return gres;
}


void run(double *TP, int *pn, double *res)
{
    *res = F(TP, pn);
    return;
}


#define PROBDIM  3
void do_mytest(double X[PROBDIM], int N, int IORD);


int main(int argc, char *argv[])
{
    int i;
    double X[PROBDIM];

    for (i = 0; i < PROBDIM; i++) X[i] = 1.0;

    c_pndl_init();
    torc_init(argc, argv, 0);

    int rank = 0;

    if (rank == 0)
    {
        int IORD = 2;
        if (argc == 2) IORD = atoi(argv[1]);

        do_mytest(X, PROBDIM, IORD);
    }


    torc_finalize();
    c_pndl_finalize();
    //    MPI_Finalize();

    return 0;
}


void do_mytest(double X[PROBDIM], int N, int IORD)
{
    double GRAD[PROBDIM];         // 1st derivatives
    double HES[PROBDIM][PROBDIM]; // hessian
    int i, j;
    double t1, t2;

    double XL[PROBDIM], XU[PROBDIM], UH[PROBDIM];
    double FEPS = 1e-6;
    int IPRINT = 0;
    int NOC;
    int IERR;

    for (i = 0; i < PROBDIM; i++) XL[i] = 0.0;
    for (i = 0; i < PROBDIM; i++) XU[i] = X[i] + 2.0;
    for (i = 0; i < PROBDIM; i++) UH[i] = 1e-3;

#if 1    // GRADIENT
    printf("\n================================================\n");
    t1 = my_gettime();
    c_pndlga(F,X,&N,XL,XU,UH,&FEPS,&IORD,&IPRINT,GRAD,&NOC,&IERR);
    t2 = my_gettime();

    if (IORD == 4) IORD = 2;
    printf("IORD = %d\n", IORD);
    printf("NOC = %d\n", NOC);
    printf("t2-t1 = %lf seconds\n", t2-t1);
    printf("GRADIENT VECTOR :\n");
    for (i = 0; i < N; i++) {
        printf("%15.8lf ", GRAD[i]);
    }
    printf("\n"); fflush(0);
#endif

#if 1    // HESSIAN WITH FUNCTION CALLS
    printf("\n================================================\n");
    t1 = my_gettime();
    c_pndlhfa(F,X,&N,XL,XU,UH,&FEPS,&IORD,&IPRINT,(double *)HES,&N,&NOC,&IERR);
    t2 = my_gettime();

    for (i = 0; i < N; i++) {
        for (j = i+1; j < N; j++)
            HES[j][i] = HES[i][j];
    }

    printf("IORD = %d\n", IORD);
    printf("NOC = %d\n", NOC);
    printf("t2-t1 = %lf seconds\n", t2-t1);
    printf("HESSIAN MATRIX :\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
            printf("%15.8lf ", HES[i][j]);
        printf("\n");
    }
    printf("\n"); fflush(0);
#endif


#if 1    // GRADIENT+HESSIAN WITH FUNCTION CALLS
    double GRAD2[PROBDIM];    // 1st derivatives
    double HES2[PROBDIM][PROBDIM];    // hessian

    printf("\n================================================\n");
    t1 = my_gettime();
    c_pndlghfa(F,X,&N,XL,XU,UH,&FEPS,&IORD,&IPRINT,GRAD2,(double *)HES2,&N,&NOC,&IERR);
    t2 = my_gettime();

    for (i = 0; i < N; i++) {
        for (j = i+1; j < N; j++)
            HES2[j][i] = HES2[i][j];
    }

    printf("IORD = %d\n", IORD);
    printf("NOC = %d\n", NOC);
    printf("t2-t1 = %lf seconds\n", t2-t1);

    printf("GRADIENT VECTOR 2:\n");
    for (i = 0; i < N; i++) {
        printf("%15.8lf ", GRAD2[i]);
    }
    printf("\n");

    printf("HESSIAN MATRIX 2:\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
            printf("%15.8lf ", HES2[i][j]);
        printf("\n");
    }
    printf("\n"); fflush(0);
#endif

}
