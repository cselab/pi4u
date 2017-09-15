# rosebrock function 

import math
import time

def fitfun0(x,n):
#	time.sleep(1)
	s = 0.0
	for i in range (1, n):
		s += 100.0*(x[i]-x[i-1]**2.0)**2.0 + (1-x[i-1])**2.0

#	print s
	return -s
