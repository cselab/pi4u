# multivariate gaussian

import numpy as np
from scipy.stats import multivariate_normal as mvnorm
import math

def fitfun1(x, dim):

#	print('x=',x)

	mean = [0, 0]
	cov = [[1, 0], [0, 1]]  # diagonal covariance

	s = mvnorm.pdf(x, mean, cov)
	s = math.log(s)

#	print('s=',s)
	return s
