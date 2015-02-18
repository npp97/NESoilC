library(MODISTools)
library(NCBI2R)

dat <- read.csv("cprj_site.csv")

names(dat)[names(dat) == "Longitude"] <- "long"
names(dat)[names(dat) == "Latitude"] <- "lat"
names(dat)[names(dat) == "Sample_end_latest"] <- "end.date"
dat$end.date <- as.Date(dat$end.date)
dat$ID <- paste("lat", sprintf("%.4f",dat$lat),"lon", sprintf("%.4f",dat$long),"End", dat$end.date, sep="")


dat_o=dat

products <- "MOD15A2"
bands <- c("FparLai_QC","Lai_1km")

if(file.exists("MOD15A2_IndividualPixel") == FALSE)
	dir.create("MOD15A2_IndividualPixel")


flst<-dir(pattern='.asc',path='.//MOD15A2_IndividualPixel')
lat<-as.numeric(substr(flst,4,10))
long<-as.numeric(substr(flst,14,21))
lat<-lat[!is.na(lat)]
long<-long[!is.na(long)]
ii<-which(!((dat$lat%in%lat)&(dat$long%in%long)))
dat<-dat[ii,]
dat1<-dat[nrow(dat):(floor(length(ii)/2-1)),]


MODISSubsets(dat1, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD15A2_IndividualPixel",TimeSeriesLength = 3)
take2 <- UpdateSubsets(dat_o, StartDate = FALSE, Dir = "MOD15A2_IndividualPixel")
MODISSubsets(take2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD15A2_IndividualPixel", TimeSeriesLength = 3)
#MODISTimeSeries(Dir = "MOD15A2_IndividualPixel", Band = "Lai_1km")
