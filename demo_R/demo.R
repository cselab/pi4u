source("../lib_R/tmcmc.R")

# main

Nth <- c(2)
MaxStages <- c(20)
PopSize <- c(4096)

lb <- c(-10, -10)
#lb <- c(  0,   0)
ub <- c( 10,  10)
id <- c(2)


Nth
MaxStages
PopSize
lb
ub

logEv <- tmcmc("fitfun2",Nth,MaxStages,PopSize,lb,ub,id)

logEv

#tmcmc_finalize()
