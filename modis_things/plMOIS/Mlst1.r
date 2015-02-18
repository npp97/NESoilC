#install.packages("devtools")
#library(devtools)
#install_github("seantuck12/MODISTools")
library(MODISTools)
require(NCBI2R)

sessionInfo()
# MODISTools_0.94.3

## setwd()

dat <- read.csv("modis_to_q.csv")

## Making column names match what is needed by MODISTools
names(dat)[names(dat) == "Longitude"] <- "long"
names(dat)[names(dat) == "Latitude"] <- "lat"
names(dat)[names(dat) == "Sample_end_latest"] <- "end.date"
dat$end.date <- as.Date(dat$end.date)
dat$ID <- paste("lat", sprintf("%.4f",dat$lat),"lon", sprintf("%.4f",dat$long),"End", dat$end.date, sep="")

dat_o=dat
## Creating my unique ID for each unique site
products <- "MOD11A2"
bands <- c("LST_Day_1km", "QC_Day","LST_Night_1km","QC_Night","Clear_sky_days","Clear_sky_nights")

if(file.exists("MOD11A2_IndividualPixel") == FALSE)
	dir.create("MOD11A2_IndividualPixel")

dir(pattern='.asc',path='.//MOD11A2_IndividualPixel')->flst
lat<-as.numeric(substr(flst,4,10))
long<-as.numeric(substr(flst,14,21))
ii<-which(!((dat$lat%in%lat)&(dat$long%in%long)))
dat<-dat[ii,]


## Call to download the subsets
iend<-floor(dim(dat)[1]/5)
il<-((1:iend)-1)*5+1
iu<-(1:iend)*5
iu[which.max(iu)]<-dim(dat)[1]
imid<-floor(iend/2)-1
imid2<-floor(iend/2)+1

for (i in 1:imid2){

	dat2<-dat[il[i]:iu[i],]	
## Call to download the subsets
	try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD11A2_IndividualPixel",TimeSeriesLength = 3))
#	k=0
	while (class(try.err) == "try-error") {
		Delay(5)
#		k=k+1
		try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD11A2_IndividualPixel", TimeSeriesLength = 3))
#		if(k>12)break;
	}
}
# To make sure that all time-series were downloaded, run the data frame through update subsets
# to create a new data frame of all the site that might have failed to download first time through
#take2 <- UpdateSubsets(dat_o, StartDate = FALSE, Dir = "MOD15A2_IndividualPixel")
# Now to use MODISSubsets again
#MODISSubsets(take2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD11A2_IndividualPixel", TimeSeriesLength = 3)
#MODISTimeSeries(Dir = "MOD15A2_IndividualPixel", Band = "Lai_1km")
#MODISTimeSeries(Dir = "MOD11A2_IndividualPixel", Band = "LST_Day_1km")
#MODISTimeSeries(Dir = "MOD11A2_IndividualPixel", Band = "LST_Night_1km")
############################################
## If all bands use the same QA band, they can be summarised together
#MODISSummaries(dat, Dir = "MOD15A2_IndividualPixel", Product = products, Bands = bands[3:4], ValidRange = c(0, 100), ScaleFactor = c(0.01,0.1), StartDate = FALSE, Interpolate = TRUE, Yield = TRUE, QualityScreen = TRUE, QualityBand = "FparLai_QC", QualityThreshold = 0, DiagnosticPlot = TRUE)

