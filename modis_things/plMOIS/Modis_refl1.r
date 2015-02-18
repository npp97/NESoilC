library(MODISTools)
library(NCBI2R)

dat <- read.csv("modis_to_q.csv")

names(dat)[names(dat) == "Longitude"] <- "long"
names(dat)[names(dat) == "Latitude"] <- "lat"
names(dat)[names(dat) == "Sample_end_latest"] <- "end.date"
dat$end.date <- as.Date(dat$end.date)
dat$ID <- paste("lat", sprintf("%.4f",dat$lat),"lon", sprintf("%.4f",dat$long),"End", dat$end.date, sep="")

dat_o=dat

products <- "MOD09A1"
bands <- GetBands('MOD09A1')


if(file.exists("MOD09A1_IndividualPixel") == FALSE)
	dir.create("MOD09A1_IndividualPixel")


flst<-dir(pattern='.asc',path='.//MOD09A1_IndividualPixel')
lat<-as.numeric(substr(flst,4,10))
long<-as.numeric(substr(flst,14,21))
ii<-which(!((dat$lat%in%lat)&(dat$long%in%long)))
dat<-dat[ii,]

iend<-floor(dim(dat)[1]/5)
il<-((1:iend)-1)*5+1
iu<-(1:iend)*5
iu[which.max(iu)]<-dim(dat)[1]
imid1<-floor(iend/4)-1
imid2<-floor(iend/2)-1
imid3<-floor(iend*0.75)-1

for (i in 1:imid1){

	dat2<-dat[il[i]:iu[i],]	
## Call to download the subsets
	try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD09A1_IndividualPixel",TimeSeriesLength = 3))
#	k=0
	while (class(try.err) == "try-error") {
		Delay(5)
#		k=k+1
		try.err<-try(MODISSubsets(dat2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD09A1_IndividualPixel", TimeSeriesLength = 3))
#		if(k>12)break;
	}
}
#take2 <- UpdateSubsets(dat_o, StartDate = FALSE, Dir = "MOD09A1_IndividualPixel")
#MODISSubsets(take2, Products = products, Bands = bands, Size = c(0,0), SaveDir = "MOD09A1_IndividualPixel", TimeSeriesLength = 3)
#MODISTimeSeries(Dir = "MOD15A2_IndividualPixel", Band = "Lai_1km")
