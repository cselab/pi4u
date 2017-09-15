# two multivariate gaussians

import numpy as np
from scipy.stats import multivariate_normal as mvnorm
import math

def fitfun2(x,n):

#	print('x=',x)

	mean1 = [-5, -5]
	cov1 = [[1, 0], [0, 1]]  # diagonal covariance

	mean2 = [5, 5]
	cov2 = [[1, 0], [0, 1]]  # diagonal covariance

	s = mvnorm.pdf(x, mean1, cov1) + mvnorm.pdf(x, mean2, cov2)
	s = math.log(s)

#	print('s=',s)
	return s
