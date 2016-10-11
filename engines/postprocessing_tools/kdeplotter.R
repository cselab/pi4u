source("kdepairs.default.R")
source("kdepairs.R")
par(cex.axis=2)

cat("Did you delete all the variables you wanted? If no, run 'rm(list=ls())'")

# ----------------------------------------------------------------------
# name of the data file
name = "curgen_db_007"

# names of variables
l = c(expression(theta[1]), expression(theta[2]))

# flag saying if we use discrepancy of likelihood; used only for colouring
discrepancy <- 0

# flag saying if we use likelihood/log-likelihood
truelik <- 0

# ----------------------------------------------------------------------

cat("Reading file", paste0(name,".txt"),"\n")
data <- read.table(paste0(name, ".txt"))
d <- dim(data)

md <- colMeans(data[1:d[2]-1])
m2d <- colMeans(data[1:d[2]-1]^2)
sd <- sqrt(m2d-md^2)
COVd <- 100 * abs(sd / md)
cat("means:", md, "\n")
cat("stds:", sd, "\n")
cat("COVs (%):", COVd, "\n")

best_id <- which.max( data[,d[2]] )
cat("params with max likelihood: \n")
print( data[best_id,] )

if(discrepancy) data[,d[2]]=-data[,d[2]]

if(truelik) data[,d[2]]=exp(data[,d[2]])

#------- this option is for big pictures ---------
# png(paste0(name, ".png"), width=1080, height=1080, units='px', res=150, pointsize=10)
# kdepairs(data, n_1d=20, n_2d=200, labels=l)
#------- this option is for small pictures -------
png(paste0(name, ".png"), width=380, height=380, units='px', res=100, pointsize=10)
kdepairs(data, n_1d=15, n_2d=200, labels=l)
#-------------------------------------------------
dev.off()
