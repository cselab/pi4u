# This function is invoked by the C program r_test
# install.packages("mvtnorm")

library(mvtnorm)

fitfun2 <- function(x) {
#        cat("R received: ", x, "\n");

        n <- length(x)

	r <- dmvnorm(x, mean=c(5,5)) + dmvnorm(x, mean=c(-5,-5))

	r <- log(r)
        return (r)
}

