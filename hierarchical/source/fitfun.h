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
#include "engine_tmcmc.h"




void fitfun_init( data_t data );

double fitfun(double *x, int n, void *output, int *info);

void fitfun_cleanup( );
