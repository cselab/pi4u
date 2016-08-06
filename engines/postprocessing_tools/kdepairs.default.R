kdepairs.default <- function(x, n=25, density=TRUE, contour=TRUE,labels, ...) {
	require(MASS)
	require(sp)
    
	ly=length(x)-1
	last=length(x)
	y <- data.frame(x[,1:ly])
	assign("loglh",data.frame(x[,last]) , envir = .GlobalEnv)
	
	fun.lower <- function(x1, x2, ...) {
		if (is.factor(x1)) {
			x1 <- as.integer(x1)
		}
		if (is.factor(x2)) {
			x1 <- as.integer(x2)
		}
		OK <- length(unique(x1))>2 && length(unique(x2))>2
		
		if (!density && !contour) n <- 0
    
		if (n>0 && OK) {
			if (density || contour) {
        d <- tryCatch({
          kde2d(x1, x2, n=n)
          }, error = function(e) {
            cat("Error caught: ")
            print(e$message)
            if (e$message == "bandwidths must be strictly positive") {
              d <- kde2d(x1, x2, n=n, h=c(0.1*(max(x1)-min(x1)),0.1*(max(x2)-min(x2))))
            } else {
              cat("Handling of this error is not implemented!")
            }
          })
			}

			if (density) {
			vv<-bpy.colors(n = 100, cutoff.tails = 0.3, alpha = 0.5)
			image(d, col=vv, add=TRUE)
			}
			
			if (contour)
			graphics:::contour(d,add=TRUE)
		} else points(x1, x2)
	}
    
	fun.upper <- function(x1, x2, ...) {
		if (is.factor(x1)) x1 <- as.integer(x1)
		if (is.factor(x2)) x1 <- as.integer(x2)

		vv<-bpy.colors(n = 100, cutoff.tails = 0.3, alpha = 0.5)
		
		# This adds a column of color values based on the y values
		length(loglh)
		loglh
		datc <- vv[as.numeric(cut(as.numeric(loglh[,1]),breaks = 100))]
		points(x1,x2,col=datc)
    }

	fun.diag.kde = function(x, ...) {
		pu <- par("usr")
		d <- density(x,...)
		par("usr" = c(pu[1:2], 0, max(d$y)*1.5))
		lines(d)
		polygon(d, col="green", border="black")
		par("usr" = c(pu[1:2], 0, 1))

		# the following 2 lines draw a box-and-whiskers plot
		# the bottom and top of the box are always the first and third quartiles,
		# and the band inside the box is always the second quartile (the median)
		# if varwidth is TRUE, the boxes are drawn with widths proportional
		# to the square-roots of the number of observations in the groups.
		# If range is positive, the whiskers extend to the most extreme data point
		# which is no more than range times the interquartile range from the box.
		h<-boxplot(x,plot=FALSE, range=1.5)
		bxp(h,at=0.1,varwidth=TRUE,boxwex=0.1,
			horizontal=TRUE,outline=FALSE,add=TRUE,axes=FALSE)
	}

	fun.diag.hist <- function(x, ...) {
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(usr[1:2], 0, 1.5) )
		h <- hist(x, breaks = 30,plot=FALSE)
		breaks <- h$breaks; nB <- length(breaks)
		y <- h$density; #y <- y/max(y)
		rect(breaks[-nB], 0, breaks[-1], y, col="green", ...)
		box()
		h<-boxplot(x, at=0.5,varwidth=TRUE, boxwex=0.3, horizontal=TRUE, outline=FALSE,plot=FALSE)
		bxp(h,at=0.5,varwidth=TRUE,boxwex=0.3,horizontal=TRUE,outline=FALSE,add=TRUE,axes=FALSE)
	}
	
	pairs.default(y,labels=labels,lower.panel=fun.lower, upper.panel=fun.upper, diag.panel=fun.diag.kde)
	invisible(NULL)
}
