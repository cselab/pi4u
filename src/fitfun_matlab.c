#include <stdio.h>
#include <pthread.h>
#include <string.h>

static char fitfun_file[256];

void fitfun_initialize(char *fitfun_name)
{
	strcpy(fitfun_file, fitfun_name);
}

void fitfun_finalize()
{
}

#include "mex.h"

double fitfun(double *x, int N, void *output, int *info)
{
     	mxArray *in[2];
        mxArray *out[1];

        in[0] = mxCreateDoubleMatrix(1,N,mxREAL);
        double *x_in = mxGetPr(in[0]);
        memcpy(x_in, x, N*sizeof(double));

        in[1] = mxCreateDoubleScalar(N);

	mexCallMATLAB(1,&out[0],2,in,fitfun_file);
//        mexCallMATLAB(1,&out[0],2,in,"fitfun");

        mxDestroyArray(in[0]);
        mxDestroyArray(in[1]);

        double *res = mxGetPr(out[0]);
        double f = *res;
        mxDestroyArray(out[0]);

        return f;
}


#include "engine_tmcmc.h"
#include "mex.h"

//extern void tmcmc(double *res, int tmcmc_info[4], int Nth, int MaxStages, int PopSize, double *lb, double *ub);

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
//	int m, n;
//	m = mxGetM(prhs[0]);
//	n = mxGetN(prhs[0]);

//	printf("m = %d\n", m);
//	printf("n = %d\n", n);


        int option = mxGetScalar(prhs[0]);

        if (option == 0)
        {
		char fitfun_name[256];

		/* Copy the string data from prhs[0] into a C string input_buf. */
		int buflen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
		/*int status =*/ mxGetString(prhs[1], fitfun_name, buflen);

		//int i = tmcmc_initialize("fitfun");
		printf("fitfun_name = %s\b", fitfun_name);
		int i = tmcmc_initialize(fitfun_name);
                plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);     /* output: res */

                double *pres = mxGetPr(plhs[0]);
                *pres = (double)i;
        }
	else if (option == 1)
        {
                tmcmc_finalize();

        }
	else if (option == 2)
        {
                double *dtmcmc_info = mxGetPr(prhs[1]);

                int i, tmcmc_info[4];
                for (i = 0; i < 4; i++) tmcmc_info[i] = dtmcmc_info[i];

                int Nth = mxGetScalar(prhs[2]);
                int MaxStages = mxGetScalar(prhs[3]);
                int PopSize = mxGetScalar(prhs[4]);
                double *lb = mxGetPr(prhs[5]);
                double *ub = mxGetPr(prhs[6]);

                plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);     /* output: res */
                double *pres = mxGetPr(plhs[0]);

                tmcmc(pres, tmcmc_info, Nth, MaxStages, PopSize, lb, ub);
        }

	return;
}

