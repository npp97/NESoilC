# TODO: Add comment
# 
# Author: lenovo
###############################################################################

require(dplR)
require(pander)
require(doBy)
require(robustbase)
# ----------------------
sumfun <- function(x, ...) {
	c(m = mean(x, ...), md = median(x, ...), v = sd(x, ...), l = length(x, ...), mi = min(x, ...), mx = max(x, ...))
}

loc.spe <- function(tmp) {
	spe <- locc <- ag <- rt<- NA
	spp <- agrep("SpeciesName", tmp)
	locp <- agrep("Location", tmp)
	agp <- agrep("Length", tmp)
	adiam <- agrep("Project",tmp)
	
	if (locp) locc <- trim.spaces(unlist(strsplit(unlist(strsplit(tmp[locp], "="))[2], "-"))[1])
	else locc<-NA;
	rt<-c(locc)

	if (spp) spe <- trim.spaces(unlist(strsplit(tmp[spp], "="))[2])
	else spe<-NA
		rt<-c(rt,spe)
	
	if (agp) ag <- as.numeric(unlist(strsplit(tmp[agp], "="))[2])
	else ag<-NA
		rt<-c(rt,ag)
		
	if (adiam){
		a= unlist(strsplit(tmp[adiam], "D"))
		diam <- as.numeric(a[2])
	}
	else diam<-NA
		rt<-c(rt,diam)

	return(rt)
}


l.nna <- function(x) {
	return(length(which(!is.na(x))))
}

to.chron.by.col <- function(x) {
	nc <- names(x)
	lvl.nc <- levels(as.factor(nc))
	l.lnc <- length(lvl.nc)
	ch.p.site <- as.data.frame(matrix(NA, ncol = l.lnc, nrow = nrow(x)))
	names(ch.p.site) <- as.character(lvl.nc)
	for (i in 1:l.lnc) {
		ii <- which(nc %in% lvl.nc[i])
		if (length(ii) < 2) 
			ch.p.site[, i] <- x[, ii] else ch.p.site[, i] <- chron(x[, ii])[1]
	}
	return(ch.p.site)
}

grw_crv <- function(pars, x, y) {
	m <- pars[1] * x^pars[1] + pars[3]
	return(sum((m - y)^2))
}

# -------------------------------------------------------------------------------------------------------------------------------------- please set the
# data fold as the working fold.  ##### Search plot age using tree-ring chronsequences ##### to determine the Age-Diam curves for each tree and each plot
# ##### by Junhui ZHANG @ 2014-05 #####

# read raw tree-ring chronsequences -------------------------------------
flst <- dir(pattern = "*.fh", recursive = TRUE, ignore.case = TRUE)
lflst <- length(flst)
# setup site_info table from tree-ring files --------------------------------------------------
las <- as.data.frame(matrix(NA, ncol = 6, nrow = lflst))
names(las) <- c("TRFN", "location", "species", "age", "diam", "filename")
# begin to read raw chronsequence and to build the cross-section of raw-chron and site-info
con <- file(flst[1], "r", blocking = FALSE)
tmp <- readLines(con, 22)
close(con)
las[1, 6] <- flst[1]
las[1, 2:5] <- loc.spe(tmp)


for (i in 2:lflst) {
	try({
				chr1.rwl <- read.rwl(flst[i], format = "auto")
				mxm <- max(as.numeric(row.names(chr1.rwl)))
				mim <- min(as.numeric(row.names(chr1.rwl)))
				if ((mxm < 2012) & (mim > 1730)) {
					# chr.rwl = combine.rwl(chr.rwl, chr1.rwl);
					con <- file(flst[i], "r", blocking = FALSE)
					tmp <- readLines(con, 22)
					close(con)
					las[i, 6] <- flst[i]
					las[i, 2:5] <- loc.spe(tmp)
					print(i)
				}
			})
}

# save(chr.rwl,file='chr_rwl.RData') plot(chron(chr.rwl))
for (i in 1:nrow(las)) {
	m <- unlist(strsplit(flst[i], "/"))
	las$TRFN[i] <- m[length(m)]
}

las$age <- as.numeric(las$age)
las$species <- as.factor(las$species)
las$TRFN <- gsub("-", "", las$TRFN)

## ---------------------------------------------------------------------# build site info by mergeing 年轮-样地对应表 and tree-ring built site_info #
# ---------------------------------------------------------------------#
#site_info <- merge(las, site_tring)
site_info <- las
l.lst <- nrow(site_info)

chr.rwl <- read.rwl(flst[1], format = "auto")
for (i in 2:l.lst) {
	try({
				chr1.rwl <- read.rwl(flst[i], format = "auto")
				mxm <- max(as.numeric(row.names(chr1.rwl)))
				mim <- min(as.numeric(row.names(chr1.rwl)))
				if ((mxm < 2012) & (mim > 1730)) {
					chr.rwl = combine.rwl(chr.rwl, chr1.rwl)
					print(i)
				}
			})
}

