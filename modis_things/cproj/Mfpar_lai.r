#install.packages("devtools")
#library(devtools)
#install_github("seantuck12/MODISTools")
library(MODISTools)
require(NCBI2R)

sessionInfo()
# MODISTools_0.94.3

## setwd()

dat <- read.csv("cprj_site.csv")

## Making column names match what is needed by MODISTools
names(dat)[names(dat) == "Longitude"] <- "long"
names(dat)[names(dat) == "Latitude"] <- "lat"
names(dat)[names(dat) == "Sample_end_latest"] <- "end.date"
dat$end.date <- as.Date(dat$end.date)
dat$ID <- paste("lat", sprintf("%.4f",dat$lat),"lon", sprintf("%.4f",dat$long),"End", dat$end.date, sep="")

dat_o=dat
## Creating my unique ID for each unique site
products <- "MOD15A2"
#bands <- c("FparExtra_QC","FparLai_QC","Lai_1km","Fpar_1km")
bands <- c("FparLai_QC","Lai_1km")

if(file.exists("MOD15A2_IndividualPixel") == FALSE)
	dir.create("MOD15A2_IndividualPixel")


dir(pattern='.asc',path='.//MOD15A2_IndividualPixel')->flst
lat<-as.numeric(substr(flst,4,11))
long<-as.numeric(substr(flst,15,23))
ii<-which(!((dat$lat%in%lat)&(dat$long%in%long)))
dat<-dat[ii,]


## Call to download the subsets
iend<-floor(nrow(dat)/5)
il<-((1:iend)-1)*5+1
iu<-(1:iend)*5
iu[which.max(iu)]<-nrow(dat)


for (i in 1:iend){

	dat2<-dat[il[i]:iu[i],]	
	print(il[i])
	print(iu[i])
	
## Call to download the subsets
	try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD15A2_IndividualPixel",TimeSeriesLength = 3))
#	k=0
	while (class(try.err) == "try-error") {
		Delay(5)
#		k=k+1
		try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD15A2_IndividualPixel", TimeSeriesLength = 3))
#		if(k>12)break;
	}
}
# To make sure that all time-series were downloaded, run the data frame through update subsets
# to create a new data frame of all the site that might have failed to download first time through
take2 <- UpdateSubsets(dat_o, StartDate = FALSE, Dir = "MOD15A2_IndividualPixel")
# Now to use MODISSubsets again
MODISSubsets(take2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD15A2_IndividualPixel", TimeSeriesLength = 3)
#MODISTimeSeries(Dir = "MOD15A2_IndividualPixel", Band = "Lai_1km")
#MODISTimeSeries(Dir = "MOD15A2_IndividualPixel", Band = "Lai_1km")
############################################
## If all bands use the same QA band, they can be summarised together
#MODISSummaries(dat, Dir = "MOD15A2_IndividualPixel", Product = products, Bands = bands[3:4], ValidRange = c(0, 100), ScaleFactor = c(0.01,0.1), StartDate = FALSE, Interpolate = TRUE, Yield = TRUE, QualityScreen = TRUE, QualityBand = "FparLai_QC", QualityThreshold = 0, DiagnosticPlot = TRUE)

