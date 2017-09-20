# data-driven inference

import numpy as np
import math

def fitfun3(x, dim):

#	print('x=',x)

	a = x[0]
	sigma = x[1]

        data = np.loadtxt('data3.txt')
        N=data.shape[0];

	fx = data[:,0]
        fy = data[:,1]
	
	y = a*fx;

	SSE = np.sum((y-fy)**2)
        sy=sigma;
	logn = -0.5*N*math.log(2*math.pi)-0.5*N*math.log(sy*sy)-0.5*SSE/(sy*sy)

#	print('sse=',SSE)
#	print('logn=',logn)

	return logn
