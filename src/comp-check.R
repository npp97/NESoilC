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


 #extract basename without extension; e.g. feed XXXX.fh and output XXXX
''

dict_sp<-read.csv('dict_sp.csv',header=F,colClasses=c('character'))
dict_site<-read.csv('dict_site.csv',header=F,colClasses=c('character'))
flst.fh<-dir(pattern='.fh',recursive=T,full.names = TRUE)
flst.rwl<-dir(pattern='.rwl',recursive=T,full.names = TRUE)
flst.all<-c(flst.fh,flst.rwl)

flst1<-basename(flst.all)
f.ext<-substr(flst1,nchar(flst1)-2,nchar(flst1))

flst<-sub("^([^.]*).*", "\\1", flst1)

site<-matrix(NA,ncol=1,nrow=length(flst))
mtt<-ydIDt<-mtd<-ydIDd<-age<-sp<-site
diam<-height<-mtt;

# tree species
for (i in 1:nrow(dict_sp)){ 
	ii<-grep(dict_sp$V1[i],substr(flst,5,nchar(flst)))
	if (length(ii)>0){sp[ii]<-dict_sp$V1[i];ii<-NULL;}
}

#mountain name 
for (j in 1:nrow(dict_site)){
	ii1<-grep(dict_site$V1[j],substr(flst,1,5))
	if (length(ii1)>0) {site[ii1]<-dict_site$V1[j];ii1<-NULL;}
}

#ydID

#length

#Extract Diam and Height from strings like filename and KeyCode/Project

ii<-grep(pattern="D\\d{1,}|D\\d{1,}.\\d{1,}",flst1) #find the filename containing Height and Diameter
regmatches(flst1[ii],regexpr("D\\d{1,}|D\\d{1,}.\\d{1,}",flst1[ii]))->D     #got the charater
sub("H\\d{1,}|H\\d{1,}.\\d{1,}","",D)->D
diam[ii]<-as.numeric(substr(D,2,nchar(D)))

ii<-grep(pattern="H\\d{1,}|H\\d{1,}.\\d{1,}",flst1) #find the filename containing Height and Diameter
regmatches(flst1[ii],regexpr("H\\d{1,}|H\\d{1,}.\\d{1,}",flst1[ii]))->H
height[ii]<-as.numeric(substr(H,2,nchar(H)))

for (i in 1:length(f.ext)){
	if(toupper(f.ext[i])=='.FH'){
		inp <- readLines(flst.all[i], ok = TRUE, warn = FALSE)
		ii<-grep(pattern="H\\d{1,}|H\\d{1,}.\\d{1,}",inp)
		if(length(ii)>0){
			H<-regmatches(inp[ii],regexpr("H\\d{1,}|H\\d{1,}.\\d{1,}",inp[ii]));
			height[i]<-as.numeric(substr(H,2,nchar(H)))
		}
		ii<-grep(pattern="D\\d{1,}|D\\d{1,}.\\d{1,}",inp)
		if(length(ii)>0){
			D<-regmatches(inp[ii],regexpr("D\\d{1,}|D\\d{1,}.\\d{1,}",inp[ii]))
			diam[i]<-as.numeric(substr(D,2,nchar(D)))
		}	
	}
}


###processing of XXAL and DXAL data
#dy<-read.csv('年轮-样地对应表.csv',colClasses=c(rep('character',18)))
#fn.full<-dir(pattern='.fh',recursive=T)
#fn.bn<-basename(fn.full)
#fn.dir<-dirname(fn.full)
#fn<-data.frame(fn.full,TRFN=fn.bn,fn.dir)
#
#dy$TRFN <- gsub("-", "", dy$TRFN)
#fn$TRFN <- gsub("-","",fn$TRFN)
#dy$TRFN <- gsub(" ", "\\1", dy$TRFN)
#fn$TRFN <- gsub(" ","\\1",fn$TRFN)
#fndy<-merge(dy,fn)
#
#new_name<-paste(as.character(fndy$fn.dir),as.character(fndy$filename.new),sep='/')
#file.rename(from=as.character(fndy$fn.full),to=new_name)
#

write.csv(cbind(flst,site,sp),file='site_sp.csv')



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


