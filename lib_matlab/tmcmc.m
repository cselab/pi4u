function logEv = tmcmc(fname,id,Nth,MaxStages,PopSize,lb,ub)

i=tmcmc_initialize(fname)

% the following is used if MPI parallel execution is used for the executable of main.m
if i>0
	return
end

t_info=[id 0 0 0];
logEv = mex_tmcmc(2, t_info,Nth,MaxStages,PopSize,lb,ub)

tmcmc_finalize()

end

	
