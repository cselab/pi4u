# This function is invoked by the C program r_test

fitfun0 <- function(x) {
#       cat("R received: ", x, "\n");

        n <- length(x)
        r <- sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)

        return (-r)
}

