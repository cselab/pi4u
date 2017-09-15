tmcmc_initialize <- function(fname) {
	out <- .C("tmcmcR_initialize",
			fname=as.character(fname))
}

tmcmc_finalize <- function() {
	out <- .C("tmcmcR_finalize")
}

tmcmc_core <- function(t_info, Nth, MaxStage, PopSize, lb, ub) {
		logEv <- 0
		out <- .C("tmcmcR",
			logEv=as.double(logEv), 
			t_info=as.integer(t_info), 
			Nth=as.integer(Nth), 
			MaxStages=as.integer(MaxStages), 
			PopSize=as.integer(PopSize), 
			lb=as.double(lb), 
			ub=as.double(ub))

		return (out$logEv)
}

tmcmc <- function(fname, Nth, MaxStage, PopSize, lb, ub, id) {
		tmcmc_initialize(fname)

		t_info <- c(id, 0, 0, 0)
		logEv <- 0
		out <- .C("tmcmcR",
			logEv=as.double(logEv), 
			t_info=as.integer(t_info), 
			Nth=as.integer(Nth), 
			MaxStages=as.integer(MaxStages), 
			PopSize=as.integer(PopSize), 
			lb=as.double(lb), 
			ub=as.double(ub))

		tmcmc_finalize()
		return (out$logEv)
}


dyn.load("../lib_R/libtmcmcR.so");
