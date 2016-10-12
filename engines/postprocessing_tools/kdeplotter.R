source("kdepairs.default.R")
source("kdepairs.R")
par(cex.axis=2)

library(tools)

cat("Did you delete all the variables you wanted? If no, run 'rm(list=ls())'")

# ----------------------------------------------------------------------
# name of the data file
fname <- "curgen_db_007.txt"

# names of variables
l = c(expression(theta[1]), expression(theta[2]))

discrepancy <- 0  # is the last column of the data discrepancy of likelihood?
truelik <- 0  # use likelihood or log-likelihood for coloring

# ----------------------------------------------------------------------

print_nicely <- function(text, x) {
  cat(text, formatC(x, digits=6, format='f', width=8), "\n", sep="\t")
}


cat("Reading file", fname,"\n")
data <- read.table(fname)
data <- array(data=unlist(data), dim=dim(data))
nd <- dim(data)[2]  # dimension of the samples

# mean
md <- colMeans(data[, 1:nd-1])
print_nicely("means:\t", md)

# standard deviation
m2d <- colMeans(data[, 1:nd-1]^2)
sd <- sqrt(m2d-md^2)
print_nicely("stds:\t", sd)

# coefficient of variation
COVd <- 100*abs(sd/md)
print_nicely("COVs (%):", COVd)

# most probable parameters
best_id <- which.max(data[, nd])
best <- data[best_id, 1:nd-1]
print_nicely("MPVs:\t", best)

if(discrepancy) data[ ,nd] <-    -data[ ,nd]
if(truelik)     data[ ,nd] <- exp(data[ ,nd])

pname <- paste0(file_path_sans_ext(fname), ".png")
#------- this option is for big pictures ---------
# png(pname, width=1080, height=1080, units='px', res=150, pointsize=10)
# kdepairs(data, n_1d=20, n_2d=200, labels=l)
#------- this option is for small pictures -------
png(pname, width=380, height=380, units='px', res=100, pointsize=10)
kdepairs(data, n_1d=15, n_2d=200, labels=l)
#-------------------------------------------------
dev.off()
