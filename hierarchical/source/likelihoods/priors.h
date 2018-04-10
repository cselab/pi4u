#ifndef _PRIORS_H_
#define _PRIORS_H_

static inline double unifpdf(double x, double l, double u);
static inline double unifpdf2(double x, double l, double len);
static inline double trnpdf(double x, double m, double s, double l, double u);
static inline double log_lognpdf(double x, double m, double s);

//inline double log_normal_pdf(double x, double mu, double sigma);


double priorHB(double *x, double *psi, int n);
double log_priorHB(double *x, double *psi, int n);

#endif
