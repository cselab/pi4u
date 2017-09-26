# data-driven inference

import numpy as np
import math

def fitfun3(theta, dim):

#	print('theta=',theta)

	a = theta[0]
	sigma = theta[1]

	data = np.loadtxt('data3.txt')
	N=data.shape[0];

	x = data[:,0]
	d = data[:,1]
	
	y = a*x;

	SSE = np.sum((y-d)**2)
	sy=sigma;
	logn = -0.5*N*math.log(2*math.pi)-0.5*N*math.log(sy*sy)-0.5*SSE/(sy*sy)

#	print('sse=',SSE)
#	print('logn=',logn)

	return logn
