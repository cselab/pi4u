#ifndef _FITFUN_H_
#define _FITFUN_H_

#include <stdio.h>
#include <stdlib.h>

#include <fcntl.h>
#include <ftw.h>
#include <math.h>

#include <unistd.h>
#include <sys/wait.h>
#include <errno.h>

#include "spawner.h"
#include "torc.h"




void fitfun_initialize();

double fitfun(double *x, int n, void *output, int *info);

void fitfun_finalize();

#endif
