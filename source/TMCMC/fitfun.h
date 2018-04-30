#ifndef _FITFUN_H_
#define _FITFUN_H_

void fitfun_initialize(char *s);

double fitfun(double *x, int n, void *output, int *info);

void fitfun_finalize();

#endif
