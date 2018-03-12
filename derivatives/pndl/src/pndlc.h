#ifndef _PNDLC_H_
#define _PNDLC_H_

void c_pndlhfa(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *hes, const int *ld, int *noc, int *ierr);
void c_pndlga(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *grad, int *noc, int *ierr);
void c_pndlghfa(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *grad, double *hes, const int *ld, int *noc, int *ierr);
void c_pndl_init(void);
void c_pndl_finalize(void);
void c_pndlhga(void (*GRD)(double *, int *, double *), double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *hes, int *ld, int *noc, int *ierr);
void c_pndlja(void (*RSD)(double *, int *, int *, double *), double *x, const int *n, const int *m, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *fj, int *ld, int *noc, int *ierr);

#endif
