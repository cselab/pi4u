import numpy as np
import ctypes
import os

def tmcmc(fitfun,dim=2,maxstages=20,popsize=1024,lowerbound=[-6,-6],upperbound=[6,6],id=0):
	# get access to the library and its subroutines

	path = os.path.dirname(os.path.realpath(__file__))
	print(path)

	#lib = ctypes.cdll.LoadLibrary('../lib_python/libtmcmc.so')
	lib = ctypes.cdll.LoadLibrary(os.path.join(path,'libtmcmc.so'))
	
	pytmcmc_initialize = lib.tmcmc_initialize
	pytmcmc_finalize = lib.tmcmc_finalize
	pytmcmc_tmcmc = lib.tmcmc

	pytmcmc_initialize(fitfun)

	# input arguments

	print('dim=',dim)
	print('maxstages=',maxstages)
	print('popsize=',popsize)

	Nth = dim
	MaxStages = maxstages
	PopSize = popsize

	t_info = np.empty(4, dtype=np.int)
	t_info[0] = id
	t_info[1] = 0
	t_info[2] = 0
	t_info[3] = 0

	lb = np.empty(Nth, dtype=np.double)
#	lb[0] = -6
#	lb[1] = -6

	for i in range(0, Nth):
		print("i=",i)
		lb[i] = lowerbound[i]


	ub = np.empty(Nth, dtype=np.double)
#	ub[0] = 6
#	ub[1] = 6

	for i in range(0, Nth):
		ub[i] = upperbound[i]

	#print 't_info: %s' % t_info
	print("lb:", lb)
	print("ub:", ub)

	#return 0
	# output argument
	logEv = np.empty(1, dtype=np.double)


	pytmcmc_tmcmc(
		ctypes.c_void_p(logEv.ctypes.data), 
		ctypes.c_void_p(t_info.ctypes.data), 
		ctypes.c_int(Nth),
		ctypes.c_int(MaxStages),
		ctypes.c_int(PopSize),
		ctypes.c_void_p(lb.ctypes.data),
		ctypes.c_void_p(ub.ctypes.data)
	)

	print("logEv:", logEv[0])

	pytmcmc_finalize()

	return logEv
