# TODO: Add comment
# 
# Author: lenovo
###############################################################################
keys<-c("KeyCode=","DateBegin=","DateEnd=","Project=","Location=","KeyNo=","SpeciesName=","Unit=")

inp <- readLines(fname, ok = TRUE, warn = FALSE)
header.begin <- grep("^HEADER:$", inp)
header.end <- grep("^DATA:(Tree|Single)$", inp)
n <- length(header.end)


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


