#ifndef _FITFUN_H_
#define _FITFUN_H_

void fitfun_initialize(char *s);
void fitfun_finalize();
double fitfun(double *x, int N, void *output, int *info);

#endif
