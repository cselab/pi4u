# This function is invoked by the C program r_test
# install.packages("mvtnorm")

#library(mvtnorm)

fitfun3 <- function(x) {
#	cat("R received: ", x, "\n");
#	n <- length(x)

	a <- x[1]
	sigma <- x[2]

	data <- read.table('data3.txt')

	fx <- data$V1
	fy <- data$V2
	N <- length(fy)

	y <- a*fx

	SSE <- sum((y-fy)**2)
        sy<-sigma
	logn <- -0.5*N*log(2*pi)-0.5*N*log(sy*sy)-0.5*SSE/(sy*sy)

#	print('sse=',SSE)
#	print('logn=',logn)

	return (logn)
}
