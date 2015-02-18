## Appendix SA1
## Script for downloading MODIS vegetation indices for analysis
## conducted in: Tuck et al. MODISTools – downloading and processing MODIS
## remotely-sensed data in R.
## Author: Helen Phillips
## Analysis done on a MAC computer (10.8.5).

library(MODISTools)
library(NCBI2R)

sessionInfo()
# MODISTools_0.94.3

## setwd()

dat <- read.csv("cprj_site.csv")

## Making column names match what is needed by MODISTools
names(dat)[names(dat) == "Longitude"] <- "long"
names(dat)[names(dat) == "Latitude"] <- "lat"
names(dat)[names(dat) == "Sample_end_latest"] <- "end.date"
dat$end.date <- as.Date(dat$end.date)
dat_o<-dat
## Creating my unique ID for each unique site
dat$ID <- paste("lat", sprintf("%.4f",dat$lat),"lon", sprintf("%.4f",dat$long),"End", dat$end.date, sep="")

dir(pattern='.asc',path='.//MOD13Q1_IndividualPixel')->flst
lat<-as.numeric(substr(flst,4,10))
long<-as.numeric(substr(flst,14,21))
lat<-lat[!is.na(lat)]
long<-long[!is.na(long)]
ii<-which(!((dat$lat%in%lat)&(dat$long%in%long)))
dat<-dat[ii,]

## Call to download the subsets
iend<-floor(nrow(dat)/5)
il<-((1:iend)-1)*5+1
iu<-(1:iend)*5
iu[which.max(iu)]<-nrow(dat)
###############################
#### EVI and NDVI DOWNLOAD ####
###############################

## Specifying my product, and associated bands
products <- "MOD13Q1"
bands <- c("250m_16_days_pixel_reliability", "250m_16_days_NDVI", "250m_16_days_EVI")

if(file.exists("MOD13Q1_IndividualPixel") == FALSE)
	dir.create("MOD13Q1_IndividualPixel")
	

for (i in iend:(floor(iend/2)-1)){

	dat2<-dat[il[i]:iu[i],]	
## Call to download the subsets
	try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD13Q1_IndividualPixel",TimeSeriesLength = 3))
#	k=0
	while (class(try.err) == "try-error") {
		Delay(10)
#		k=k+1
		try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD13Q1_IndividualPixel", TimeSeriesLength = 3))
#		if(k>12)break;
	}
}
# T
# To make sure that all time-series were downloaded, run the data frame through update subsets
# to create a new data frame of all the site that might have failed to download first time through
#take2 <- UpdateSubsets(dat_o,StartDate = FALSE, Dir = "MOD13Q1_IndividualPixel")


# Now to use MODISSubsets again
#MODISSubsets(take2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD13Q1_IndividualPixel", TimeSeriesLength = 3)

#MODISTimeSeries(Dir = "MOD13Q1_IndividualPixel", Band = "250m_16_days_EVI")


############################################
## If all bands use the same QA band, they can be summarised together
#MODISSummaries(dat, Dir = "MOD13Q1_IndividualPixel", Product = products, Bands = bands[2:3], ValidRange = c(-2000, 10000), NoDataFill = -3000, ScaleFactor = 0.0001, StartDate = FALSE, Interpolate = TRUE, Yield = TRUE, QualityScreen = TRUE, QualityBand = "250m_16_days_pixel_reliability", QualityThreshold = 0, DiagnosticPlot = TRUE)
