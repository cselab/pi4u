# This function is invoked by the C program r_test
# install.packages("ggplot2")

library(ggplot2)
library(RColorBrewer)
#library(ggthemes)

display_gen <- function(datafile,Nth,i,j) {


data <- read.table(datafile)
dd <- data
data <- as.matrix(data)

#dd$V1 <- data$V1
#dd$V2 <- data$V2
#dd$V3 <- data$V3

dd$V1 <- data[,i]
dd$V2 <- data[,j]
dd$V3 <- data[,Nth+1]

#z<-data[,Nth+1]
#normalized = (z-min(z))/(max(z)-min(z))
#to<-c(1)
#from<-c(-1)
#normalized = (z - min(z)) / max(z - min(z)) * (to - from) + from
#dd$V3 <- normalized

#ggplot(dd, aes(x=V1,y=V2,color=V3)) + geom_point() + scale_color_gradient2(low = 'red', high = 'blue')

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min(dd$V3), max(dd$V3)))
ggplot(dd, aes(x=V1,y=V2,color=V3)) + geom_point(size=1) + sc

}
