#ifndef _PNDLC_H_
#define _PNDLC_H_

#include "pndl_config.h"

void F77_FUNC(pndlhfa,PNDLHFA)(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *hes, const int *ld, int *noc, int *ierr);
void F77_FUNC(pndlga,PNDLGA)(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *grad, int *noc, int *ierr);
void F77_FUNC(pndlghfa,PNDLGHFA)(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *grad, double *hes, const int *ld, int *noc, int *ierr);
void F77_FUNC(pndl_init,PNDL_INIT)(void);
void F77_FUNC(pndl_finalize,PNDL_FINALIZE)(void);
void F77_FUNC(pndlhga,PNDLHGA)(void (*GRD)(double *, int *, double *), double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *hes, int *ld, int *noc, int *ierr);
void F77_FUNC(pndlja,PNDLJA)(void (*RSD)(double *, int *, int *, double *), double *x, const int *n, const int *m, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *fj, int *ld, int *noc, int *ierr);

void c_pndlhfa(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *hes, const int *ld, int *noc, int *ierr)
{
    F77_FUNC(pndlhfa,PNDLHFA)(F, x, n, xl, xu, uh, feps, iord, iprint, hes, ld, noc, ierr);
}

void c_pndlga(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *grad, int *noc, int *ierr)
{
    F77_FUNC(pndlga,PNDLGA)(F, x, n, xl, xu, uh, feps, iord, iprint, grad, noc, ierr);
}

void c_pndlghfa(double (*F)(double *, int *), const double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *grad, double *hes, const int *ld, int *noc, int *ierr)
{
    F77_FUNC(pndlghfa,PNDLGHFA)(F, x, n, xl, xu, uh, feps, iord, iprint, grad, hes, ld, noc, ierr);
}

void c_pndl_init(void)
{
    F77_FUNC(pndl_init,PNDL_INIT)();
}

void c_pndl_finalize(void)
{
    F77_FUNC(pndl_finalize,PNDL_FINALIZE)();
}

void c_pndlhga(void (*GRD)(double *, int *, double *), double *x, const int *n, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *hes, int *ld, int *noc, int *ierr)
{
    F77_FUNC(pndlhga,PNDLHGA)(GRD, x, n, xl, xu, uh, feps, iord, iprint, hes, ld, noc, ierr);
}

void c_pndlja(void (*RSD)(double *, int *, int *, double *), double *x, const int *n, const int *m, const double *xl, const double *xu, const double *uh, const double *feps, const int *iord, const int *iprint, double *fj, int *ld, int *noc, int *ierr)
{
    F77_FUNC(pndlja,PNDLJA)(RSD, x, n, m, xl, xu, uh, feps, iord, iprint, fj, ld, noc, ierr);
}

#endif
