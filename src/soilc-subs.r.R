# TODO: Add comment
# 
# Author: lenovo
###############################################################################


function(ddata, bins, resolution){
#
# Function to create a histogram of non-zero data for use in 
# the mix.dist program
#
# !! REQUIRES THAT THE PACKAGE "mixdist" IS ATTACHED !!
# 
# Call: mix.dist.hist <- make.mix.data(ddata, bars, resolution)
#
# ddata: array of single value, non-zero, variables
# bins: Number of groups in the histogram
#
# data.out: A two dimensional array of histogram-type grouped data
# with first column the upper bound value of the interval and
# second column the number of values in the interval
#
	data.out <- matrix(0, bins, 2)
	nvals <- length(ddata)
#-----------------------
# Create histogram intervals
#
	max.ddata <- max(ddata)
	min.ddata <- min(ddata)
	offset <- (max.ddata - min.ddata)/resolution
	a <- min.ddata
	inter <- ((offset + max.ddata) - (min.ddata-offset))/bins
	for (i in 1:bins){
		a <- a + inter
		data.out[i, 1] <- a
	}
#
# Partition the data into bins
#
	for (i in 1:nvals){
		k <- 1
		calc.x <- ddata[i] - min.ddata
		while (calc.x > inter){
			calc.x <- calc.x - inter
			k <- k + 1
		}
		data.out[k, 2] <- data.out[k, 2] +1
	}
	data.out <- as.mixdata(as.data.frame(data.out))
	data.out 
}

sumfun <- function(x, ...){
	c(m=mean(x, ...), md=median(x,...),v=sd(x, ...), l=length(x))
}

fst2maxi <- function(x, ...){
	c(which.max(x,...),which.max(x[-which.max(x)]))
}

fst2max <- function(x, ...){
	c(max(x),max(x[-which.max(x)]))
}

