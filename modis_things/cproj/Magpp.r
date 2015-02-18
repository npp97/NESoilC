# To add two features: try-error &  lat/long check
library(MODISTools)
require(NCBI2R)  #For Delay

sessionInfo()
# MODISTools_0.94.3

## setwd()

dat<- read.csv("cprj_site.csv")

if(file.exists("MOD17A3_IndividualPixel") == FALSE)
	dir.create("MOD17A3_IndividualPixel")


## Making column names match what is needed by MODISTools
names(dat)[names(dat) == "Longitude"] <- "long"
names(dat)[names(dat) == "Latitude"] <- "lat"
names(dat)[names(dat) == "Sample_end_latest"] <- "end.date"
dat$end.date<-gsub("4","1",dat$end.date)
dat$end.date <- as.Date(dat$end.date)
dat$ID <- paste("lat", sprintf("%.4f",dat$lat),"lon", sprintf("%.4f",dat$long),"End", dat$end.date, sep="")

dat_o<-dat
products <- "MOD17A3"
bands <- c("Gpp_1km","Npp_1km","Gpp_Npp_QC_1km")

dir(pattern='.asc',path='.//MOD17A3_IndividualPixel')->flst
lat<-as.numeric(substr(flst,4,10))
long<-as.numeric(substr(flst,14,21))
ii<-which(!((dat$lat%in%lat)&(dat$long%in%long)))
dat<-dat[ii,]

l.seg<-20
## Call to download the subsets
iend<-floor(dim(dat)[1]/l.seg)
il<-((1:iend)-1)*l.seg+1
iu<-(1:iend)*l.seg
iu[which.max(iu)]<-dim(dat)[1]

for (i in 1:iend){

	dat2<-dat[il[i]:iu[i],]
	try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD17A3_IndividualPixel", TimeSeriesLength = 3))
	k=0
	while (class(try.err) == "try-error") {
#		take2 <- UpdateSubsets(dat2,StartDate = FALSE, Dir = "MOD17A3_IndividualPixel")
# Now to use MODISSubsets again
#	    try.err<-try(MODISSubsets(take2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD17A3_IndividualPixel", TimeSeriesLength = 5))

		Delay(5)
		try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD17A3_IndividualPixel", TimeSeriesLength = 3))
#		k=k+1
#		if(k>12)break;
	}

}
# To make sure that all time-series were downloaded, run the data frame through update subsets
# to create a new data frame of all the site that might have failed to download first time through
#
############################################
## If all bands use the same QA band, they can be summarised together
take2 <- UpdateSubsets(dat_o,StartDate = FALSE, Dir = "MOD17A3_IndividualPixel")
# Now to use MODISSubsets again
MODISSubsets(take2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD17A3_IndividualPixel", TimeSeriesLength = 3)
#MODISTimeSeries(Dir = "MOD17A3_IndividualPixel", Band = "Npp_1km")
#MODISSummaries(dat2, Dir = "MOD17A3_IndividualPixel", Product = products, Bands = bands[1:2], ValidRange = c(0, 65500), NoDataFill = -9, ScaleFactor = 0.0001, StartDate = FALSE, Interpolate = TRUE, Yie