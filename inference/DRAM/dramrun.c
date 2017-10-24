#include "gsl_headers.h"
#include "auxil.c"
#include <math.h>
#ifndef min
#define min(a,b) (a)<(b)?(a):(b)
#endif

#if 1
//%COVUPD covariance update
//% [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)
//function [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)
//% Marko Laine <Marko.Laine@Helsinki.FI>
//%COVUPD covariance update
//% [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)
//function [xcov,xmean,wsum]=covupd(x,w,oldcov,oldmean,oldwsum)
//% Marko Laine <Marko.Laine@Helsinki.FI>
void covupd(double *x, int start, int end, int npar, double w, double *xcov, double *xmean, double *wsum)
{
	int n = end - start;
	int p = npar;

	if (n == 0)	//% nothing to update with
	{
		//xcov = oldcov; xmean = oldmean; wsum = oldwsum;
		return;
	}
	
	//if nargin<2 | isempty(w)
	//	w = 1;
	//end
	//if length(w) == 1
	//	w = ones(n,1)*w;
	//end


	if (start > 0)	//	~isempty(oldcov)
	{ 
		//if nargin>2 & ~isempty(oldcov) % update

		double oldcov[p*p], oldmean[p], oldwsum;

		memcpy(oldcov, xcov, p*p*sizeof(double));
		memcpy(oldmean, xmean, p*sizeof(double));
		oldwsum = *wsum;


		for (int i=0; i<n; i++)
		{
			double xi[p];
			double xmeann[p];

			memcpy(xi, &x[i*p], p*sizeof(double));	//xi = x(i,:);
			*wsum   = w;

			memcpy(xmeann, xi, p*sizeof(double));	// xmeann = xi;
			//xmean  = oldmean + (*wsum)/(*wsum+oldwsum)*(xmeann-oldmean);
			for (int j=0; j<p; j++)
			{
				xmean[j] = oldmean[j] + (*wsum)/(*wsum+oldwsum)*(xmeann[j]-oldmean[j]);
			}

//			xcov =  oldcov + (*wsum)./(*wsum+oldwsum-1) .* (oldwsum/(*wsum+oldwsum) .* ((xi-oldmean)' *(xi-oldmean))  - oldcov);
			double fac = (*wsum)/(*wsum+oldwsum-1) * (oldwsum/(*wsum+oldwsum));
			double minus[p];
			double prod[p*p];
			for (int j=0; j<p; j++) minus[j]=xi[j]-oldmean[j];

			for (int j=0; j<p; j++) 
			for (int k=0; k<p; k++) 
				prod[j*p+k]=minus[j]*minus[k];


			for (int j=0; j<p; j++) 
			for (int k=0; k<p; k++) 
				xcov[j*p+k]=oldcov[j*p+k] + fac * (prod[j*p+k]-oldcov[j*p+k]);
    
			*wsum    = *wsum+oldwsum;
			memcpy(oldcov, xcov, p*p*sizeof(double));
			memcpy(oldmean, xmean, p*sizeof(double));
			oldwsum = *wsum;
		}
	}
	else
	{
		//% no update

		*wsum  = w*n; //sum(w);
		//xmean = zeros(1,p);
		//xcov  = zeros(p,p);
		memset(xmean, 0, p*sizeof(double));
		memset(xcov, 0, p*p*sizeof(double));

		for (int j=0; j<p; j++)
		{
			//xmean[j] = sum(x(:,i).*w)./(*wsum);

			double s = 0;
			for (int i = 0; i < n; i++) 
				s += x[i*npar+j]*w;

			xmean[j] = s/(*wsum);
		}

		if (*wsum>1)
		{
			//%%% (wsum-oldwsum/wsum)
			for (int i=0; i<p; i++)
			{
				for (int j=0; j<=i; j++)
				{
//					xcov(i,j) = (x(:,i)-xmean(i))' * ((x(:,j)-xmean(j)).*w)./(*wsum-1);
					double s = 0;
					for (int k = 0; k<n; k++)
					{
						s += (x[k*p+i]-xmean[i])*(x[k*p+j]-xmean[j])*w;
					}
					xcov[i*p+j] = s/(*wsum-1);

					if (i != j)
					{
						xcov[j*p+i] = xcov[i*p+j];
					}
				}
			}
		}
	}
}
#endif

void inv(double *Ainv, double *A, int n)
{
	int s;
	gsl_matrix *work = gsl_matrix_alloc(n,n), *winv = gsl_matrix_alloc(n,n);
	gsl_permutation *p = gsl_permutation_alloc(n);

	memcpy( work->data, A, n*n*sizeof(double) );
#if DEBUG
	gsl_matrix_fprintf(stdout, work, "%f");
#endif
	gsl_linalg_LU_decomp( work, p, &s );
	gsl_linalg_LU_invert( work, p, winv );
#if DEBUG
	gsl_matrix_fprintf(stdout, winv, "%f");
#endif
	memcpy( Ainv, winv->data, n*n*sizeof(double) );
	gsl_matrix_free( work );
	gsl_matrix_free( winv );
	gsl_permutation_free( p );
}

int chol(double *Achol, double *A, int n)
{
	gsl_matrix *work = gsl_matrix_alloc(n,n);
	memcpy(work->data, A, n*n*sizeof(double) );

	gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
	int status = gsl_linalg_cholesky_decomp(work);
	if (status == GSL_EDOM) {
		gsl_matrix_free(work);
		gsl_set_error_handler (old_handler); 
		return 1;	// not positive-definite (singular)
	}

#if DEBUG
	gsl_matrix_fprintf(stdout, work, "%f");
#endif
	for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		if (j<i)
			gsl_matrix_set(work, i, j, 0);

	memcpy( Achol, work->data, n*n*sizeof(double) );

#if DEBUG
	gsl_matrix_fprintf(stdout, work, "%f");
#endif

	gsl_set_error_handler (old_handler); 
	gsl_matrix_free(work);

	return 0;	// all good
}


int check_bounds(double *x, double *lbounds, double *ubounds, int n)
{
	for (int i = 0; i < n; i++)
	{
		if ((x[i]<lbounds[i]) || (x[i]>ubounds[i]))
			return 1;
	}
	return 0;
}


//function [results,chain,s2chain]=dramrun(model,data,params,options)
void dramrun(int npar)
{
//DRAMRUN  Metropolis-Hastings MCMC run with adaptive delayed rejection (DRAM)
//
// This function generates MCMC chain using DRAM adaptation for a model defined
// by user supplied sum-of-squares function and with additive i.i.d. Gaussian
// errors for the observations. The error variance sigma2 is updated
// using conjugate inverse gamma distribution.
//
// [results,chain,s2chain]=dramrun(model,data,params,options)
//
// input:
//
// model.ssfun =    ; // sum-of-squares function, ss=ssfun(par,data),
//                    // that returns  -2*log(p(y|par))
// model.priorfun = ; // prior "sum-of-squares", priorfun(par,params),
//                    // that returns -2*log(p(par)),
//                    // default: inline('0','x','params')
//
// data   = ;         // extra argument for ssfun (to pass the data etc.)
//
// params.par0   =  ; // initial parameter vector (a row vector)
// params.sigma2 =  1;// initial/prior value for the Gaussian error variance
// params.n0     = -1;// precision of sigma2 as imaginative observations
//                    //   if n0<0, no sigma2 update
// params.n      = ;  // number of actual observations (for sigma2 update)
// params.bounds = ;  // 2*npar matrix of parameter bounds
//                    // default: [-Inf,Inf]
//
// options.nsimu  = 2000;   // length of the chain
// options.qcov   = ;       // proposal covariance matrix
//
// parameters for DRAM
// options.adaptint = 10;  // how often to adapt, if zero, no adaptation
// options.drscale  = 3;   // scale for the second proposal, if zero, no DR
//
// output:
//
// results  structure that contains some info about the run
// chain    nsimu*npar MCMC chain
// s2chain  sigma² chain (if generated)


// calls covupd.m for covariance update and (optionally) gammar_mt.m for
// gamma variates

// this is a 'simple' version for demonstration and educational purposes

// Marko Laine <Marko.Laine@Helsinki.FI>
// $Revision: 1.0 $  $Date: $

	//// get values from the input structs
	int nsimu  = options.nsimu;	//getpar(options,'nsimu',10000);

	// initial parameter vector
	//par0   = getpar(params,'par0'); par0=par0(:)'; // row vector
	double *par0 = params.par0;


	// number of parameters
	//npar   = length(par0);


	// 2*npar matrix of parameter bounds
	//bounds = getpar(params,'bounds',(ones(npar,2)*diag([-Inf,Inf]))');
	double *lbounds = params.lbounds;
	double *ubounds = params.ubounds;

	if ((lbounds == NULL)||(ubounds == NULL))
	{
		printf("bounds not initialized\n");
		exit(1);
	}

	// peh: defined in normaltest.c
	// sum-of-squares function, ssfun(par,data),  -2*log(p(y|theta))
	//ssfun  = getpar(model,'ssfun');
	// prior "sum-of-squares", -2*log(p(theta))
	//
	//priorfun = getpar(model,'priorfun',inline('0','x','params'));


	////// parameters for DRAM
	// how often to adapt, if zero, no adaptation
	int adaptint = options.adaptint; 	//getpar(options,'adaptint',100);

	// scale for the second proposal, if zero, no DR
	double drscale  = options.drscale; //getpar(options,'drscale',3);

	// scale for adapting the propsal
	double adascale = options.adascale; 	//getpar(options,'adascale',2.4/sqrt(npar));

	// blow factor for covariace update
	double qcovadj  = options.qcovadj;	//getpar(options,'qcovadj',1e-5);


	// precision of sigma2 as imaginative observations
	//  if n0<0, no sigma2 update
	int n0  = params.n0;	//getpar(params,'n0',-1);
	// initial/prior value for the Gaussian error variance
	double sigma2 = params.sigma2;	//getpar(params,'sigma2',1);

#if 0
	// number of observations (needed for sigma2 update)
	int n = 0;
	if (n0 >= 0)
	{
		n = params.n;
	}
#endif

	double *qcov = options.qcov;	//getpar(options,'qcov'); // proposal covariance

	// to DR or not to DR
	int dodr;
	if (drscale<=0)
		dodr=0; 
	else
		dodr=1;

//	dodr = 0;

	int printint  = options.printint;	//getpar(options,'printint',500);
	int verbosity = options.verbosity;	//getpar(options,'verbosity',1);

	double R[npar*npar];
	double R2[npar*npar];
	double iR[npar*npar];

	chol(R, qcov, npar); // R = chol(qcov); // *adascale; // Cholesky factor of proposal covariance
	if (dodr)
	{
		for (int i = 0; i < npar*npar; i++) R2[i] = R[i]/drscale;	//R2= R./drscale; // second proposal for DR try
		inv(iR, R, npar);
	}

	double *chain = calloc(1, nsimu*npar*sizeof(double));	//  = zeros(nsimu,npar);  // we store the chain here

	double s20 = 0;
	double *s2chain;
	if (n0>=0)
	{
		s2chain = calloc(1, nsimu*sizeof(double));   // the sigma2 chain
		s20 = sigma2;
	}
	else
	{
		s2chain = NULL;
	}

	//oldpar       = par0(:)';                // first row of the chain
	double oldpar[npar];
	memcpy(oldpar, par0, npar*sizeof(double));

	//oldss        = feval(ssfun,oldpar,data);// first sum-of-squares
	//oldprior     = feval(priorfun,oldpar,params);
	double oldss    = ssfun(oldpar, npar);
	double oldprior = priorfun(oldpar, npar);

#if DEBUG
	printf("oldss = %f\n", oldss);
#endif

	int acce         = 1;                       //  how many accepted moves

	//chain(1,:)   = oldpar;
	memcpy(&chain[0*npar], oldpar, npar*sizeof(double)); 

	if (s20>0)
	{
		//s2chain(1,:) = sigma2;
		s2chain[0] = sigma2;
	}


#if 1
	// covariance update uses these to store previous values
	double chaincov[npar*npar];
	double chainmean[npar];
	double wsum = 0;
	int lasti = 0;
	////// the simulation loop
	int isimu;
	for (isimu=1; isimu<nsimu; isimu++)	// 2 .. nsimu
	{
//		if ((isimu+1)/printint == (int)((isimu+1)/printint)) // info on every printint iteration
		if ((isimu+1)%printint == 0) 
		{
			printf("isimu=%d, %d%% done, accepted: %d%%\n", isimu,(int)((float)(isimu+1)/nsimu*100),(int)(((float)acce/isimu)*100));
		}
  
		double newpar[npar];

		//newpar = oldpar+randn(1,npar)*R;     // a new proposal
		for (int i = 0; i < npar; i++)
		{
			newpar[i] = oldpar[i];
			for (int j = 0; j < npar; j++)
				newpar[i] += normalrand(0,1)*R[i*npar+j]; 
		}

		int accept;
		accept = 0;
		// check bounds
		double newss, newprior, alpha12;
		if (check_bounds(newpar, lbounds, ubounds, npar))
		{
			newss = DBL_MAX; //+Inf;
			newprior = 0;
			alpha12 = 0;
		}
		else // inside bounds, check if accepted
		{
			newss  = ssfun(newpar,npar);   // sum-of-squares
			newprior = priorfun(newpar,npar); // prior ss
			alpha12 = min(1,exp(-0.5*(newss-oldss)/sigma2-0.5*(newprior-oldprior)));
			double r = uniformrand(0,1);
			if (r < alpha12) // we accept
			{
				accept   = 1;
				acce     = acce+1;
#if DEBUG
				printf("accepting (%f<%f): from [%f %f]=(%f) to [%f %f]=(%f)\n", r, alpha12, oldpar[0], oldpar[1], oldss, newpar[0], newpar[1], newss); 
#endif
				memcpy(oldpar, newpar, npar*sizeof(double));  // oldpar = newpar;
				oldss    = newss;
				oldprior = newprior;
			}
			else
			{
#if DEBUG
				printf("not accepting (%f>%f): from [%f %f]=(%f) to [%f %f]=(%f)\n", r, alpha12, oldpar[0], oldpar[1], oldss, newpar[0], newpar[1], newss); 
#endif
			}
		}

		if (accept == 0 && dodr) // we reject, but make a new try (DR)
		{
			//newpar2 = oldpar+randn(1,npar)*R2;  // a new try
			double newpar2[npar];
			for (int i = 0; i < npar; i++)
			{
				newpar2[i] = oldpar[i];
				for (int j = 0; j < npar; j++)
					newpar2[i] += normalrand(0,1)*R2[i*npar+j]; 
			}

			double newss2, newprior2;
			if (check_bounds(newpar2, lbounds, ubounds, npar))
			{
				newss2 = DBL_MAX; //+Inf;
				newprior2 = 0;
			}
			else // inside bounds
			{
				newss2    = ssfun(newpar2,npar);
				newprior2 = priorfun(newpar2,npar);
				double alpha32 = min(1,exp(-0.5*(newss-newss2)/sigma2 -0.5*(newprior-newprior2)));
				double l2 = exp(-0.5*(newss2-oldss)/sigma2 - 0.5*(newprior2-oldprior));
				//double q1 = exp(-0.5*(norm((newpar2-newpar)*iR)^2-norm((oldpar-newpar)*iR)^2));
				double q1;
				{
					double term1[npar], minus1[npar];
					for (int i = 0; i < npar; i++) minus1[i] = newpar2[i]-newpar[i];
					for (int i = 0; i < npar; i++) 
					{
						for (int j = 0; j < npar; j++)
							term1[i] = minus1[i]*iR[i*npar+j];
					}

					double term2[npar], minus2[npar];
					for (int i = 0; i < npar; i++) minus2[i] = oldpar[i]-newpar[i];
					for (int i = 0; i < npar; i++) 
					{
						for (int j = 0; j < npar; j++)
							term2[i] = minus2[i]*iR[i*npar+j];
					}
					double norm1 = norm(term1, npar);
					double norm2 = norm(term2, npar);
					q1 = exp(-0.5*(norm1*norm1-norm2*norm2));
				}

				double alpha13 = l2*q1*(1-alpha32)/(1-alpha12);
				if (uniformrand(0,1) < alpha13) // we accept
				{
					accept = 1;
					acce     = acce+1;
					memcpy(oldpar, newpar2, npar*sizeof(double));  // oldpar = newpar2;
					oldss    = newss2;
					oldprior = newprior2;
				}
			}
 		}
  
		memcpy(&chain[isimu*npar], oldpar, npar*sizeof(double)); 	//chain(isimu,:) = oldpar; 
		// update the error variance sigma2
		if (s20 > 0)
		{
			//sigma2  = 1./gammar_mt(1,1,(n0+n)./2,2./(n0*s20+oldss));
			//s2chain(isimu,:) = sigma2;
		}
  
		if (adaptint>0 && ((isimu+1)%adaptint == 0))
		{
			// adapt the proposal covariances
			if (verbosity) printf("adapting\n");

#if 1
			// update covariance and mean of the chain
			//[chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,:),1,chaincov,chainmean,wsum);
			covupd(&chain[lasti*npar], lasti, isimu, npar, 1, chaincov, chainmean, &wsum);
			lasti = isimu;

			double chaincov_tmp[npar*npar];
			for (int i = 0; i < npar; i++)
			for (int j = 0; j < npar; j++)
				if (i == j)
					chaincov_tmp[i*npar+j] = chaincov[i*npar+j] + qcovadj;
				else
					chaincov_tmp[i*npar+j] = chaincov[i*npar+j];

			
			//[Ra,is] = chol(chaincov + eye(npar)*qcovadj);
			double Ra[npar*npar];
			int is = chol(Ra, chaincov_tmp, npar);

			if (is) // singular cmat
			{
				printf("Warning cmat singular, not adapting\n");
			}
			else
			{
				for (int i = 0; i < npar; i++)
				for (int j = 0; j < npar; j++)
					R[i*npar+j] = Ra[i*npar+j]*adascale;

				if (dodr)
				{  
					for (int i = 0; i < npar*npar; i++) R2[i] = R[i]/drscale;	//R2= R./drscale; // second proposal for DR try
					inv(iR, R, npar);
				}
			}
#endif
		}
	}
#endif


//	printf("chain:\n");
//	for (int isimu=0; isimu<nsimu; isimu++)
//	{
//		printf("%f %f\n", chain[isimu*npar+0], chain[isimu*npar+1]);
//	}

	printf("acceptance = %d\n", (int)((100.0*acce)/nsimu));
	FILE *fp = fopen(params.filename, "w");
	for (int isimu=0; isimu<nsimu; isimu++)
	{
		fprintf(fp, "%f %f\n", chain[isimu*npar+0], chain[isimu*npar+1]);
	}
	fclose(fp);

	printf("lasti: %d\n", lasti);
	printf("isimu: %d\n", isimu);

	covupd(&chain[lasti*npar], lasti, isimu, npar, 1, chaincov, chainmean, &wsum);

	for (int i = 0; i < npar; i++)
		printf("chainmean[%d]=%.4f\n", i, chainmean[i]);

	for (int i = 0; i < npar; i++)
	for (int j = 0; j < npar; j++)
		printf("chaincov[%d,%d]=%.4f\n", i, j, chaincov[i*npar+j]);

	printf("wsum = %f\n", wsum);
#if 0
	// calculate covariance and mean of the chain
	[chaincov,chainmean,wsum] = covupd(chain((lasti+1):isimu,:),1,chaincov,chainmean,wsum);
#endif

#if 0
	results.class = 'MCMC results';
	results.accepted=acce./nsimu;              // acceptance ratio
	results.mean = chainmean;
	results.cov  = chaincov;
//	results.qcov = R'*R;
	results.R = R;
	results.nsimu = nsimu;
	results.drscale = drscale;
	results.adascale = adascale;
	results.adaptint = adaptint;
#endif

}

