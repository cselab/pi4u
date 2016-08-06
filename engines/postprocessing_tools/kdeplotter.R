# ----------------------------------------------------------------------
# path to the helper files
source("~/workspace/source_codes/UQ/plot_distribution/kdepairs.default.R")
source("~/workspace/source_codes/UQ/plot_distribution/kdepairs.R")
par(cex.axis=2)

cat("Did you delete all the variables you wanted to delete? If no, run 'rm(list=ls())'")

# ----------------------------------------------------------------------
# name of the data file witout extension
name = "~/Desktop/curgen_db_000"

# names of the variables
l = c(expression(epsilon),expression(sigma),expression(p),expression(sigma[n]))

# flag saying if we use discrepancy of likelihood; used only for colouring
discrepancy <- 0
# flag saying if we use likelihood/log-likelihood
truelik <- 0

# ----------------------------------------------------------------------

cat("Reading file", paste0(name,".txt"),"\n")
data<-read.table(paste0(name, ".txt"))
d <- dim(data)

md <- colMeans(data[1:d[2]-1])
m2d <- colMeans(data[1:d[2]-1]^2)
sd <- sqrt(m2d-md^2)
COVd <- 100 * sd / md
cat("means:", md, "\n")
cat("stds:", sd, "\n")
cat("COVs (%):", COVd, "\n")

best_id <- which.max( data[,d[2]] )
cat("params with max likelihood: \n")
print( data[best_id,] )

if(discrepancy) data[,d[2]]=-data[,d[2]]

if(truelik) data[,d[2]]=exp(data[,d[2]])

# plotting device
png(paste0(name, ".png"))
#pdf(paste0(name, ".pdf"))
#tikz(paste0(name, ".tikz"))
#cairo_ps(file=paste0(name,".eps"))

kdepairs(data,n=24,labels=l)

dev.off()

# ----------------------------------------------------------------------
