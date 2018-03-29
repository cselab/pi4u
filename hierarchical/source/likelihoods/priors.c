#include "priors.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <math.h>


//TODO everything in log scale TODO



inline double unifpdf(double x, double l, double u) {
    if ((l <= x) && (x <= u))
        return 1.0/(u-l);
    return 0;
}



inline double unifpdf2(double x, double l, double len) {
    if ((l <= x) && (x <= l+len))
        return 1.0/len;
    return 0;
}



inline double lognpdf(double x, double m, double s) {
    if (s > 0)
        return gsl_ran_lognormal_pdf(x, m, s);
    else
        return NAN;
}


inline double trnpdf(double x, double m, double s, double l, double u) {
    if (s > 0)
        return gsl_ran_gaussian_pdf(fabs(x-m), s) /
            (1.0 - gsl_cdf_gaussian_P(fabs(m-l), s)
                 - gsl_cdf_gaussian_P(fabs(m-u), s));
    else
        return NAN;
}



inline double log_normal_pdf(double x, double mu, double sigma){
	if(sigma>0){
		double tmp = (x-mu)/sigma;
		return  -0.5*log(2*M_PI*sigma*sigma) - 0.5*tmp*tmp ;
	}
	else
		return NAN;
}





double log_priorHB(double *x, double *psi, int n) {

    double res = 0.0;

    for (int i = 0; i < n; i++) {
        res += log_normal_pdf( x[i], psi[2*i], psi[2*i+1] );
    }

	// printf("%lf  |   %lf  %lf  | %lf \n",x[i], psi[2*i], psi[2*i+1],   log_normal_pdf( x[i], psi[2*i], psi[2*i+1] ) );
	// printf("============= \n");

    return res;
}
