# This function is invoked by the C program r_test
# install.packages("mvtnorm")

library(mvtnorm)

fitfun1 <- function(x) {
#        cat("R received: ", x, "\n");

        n <- length(x)

	r <- dmvnorm(x, mean=c(0,0))

	r <- log(r)
        return (r)
}

