require(mgcv)
require(extremevalues)

#Beautiful usage style: rv_outliers(rv_outliers(rv_outliers(b[,4],b[,1],0.05),b[,1],0.05),b[,1],0.05)
#rv_outliers(rv_outliers(rv_outliers(b[,4],b[,1],0.05),b[,1],0.05),b[,1],0.05)->cc;points(b[,4],pch=14,col=2);abline(h=quantile(cc,c(0.05,0.95)))
#quantile(cc,c(0.05,0.95))
#rv_outliers(a[,4],b[,1],0.05)->cc;points(a[,4],pch=14,col=2);abline(h=quantile(cc,c(0.05,0.95),na.rm=T));abline(h=quantile(a[,4],c(0.05,0.95),na.rm=T),col=2)

rv_outliers<-function(y2ana,ts,tol=0.05){
	dta<-data.frame(y2ana,ts);
	ii_av<-which(!is.na(dta$y2ana));
	ii_na<-which(is.na(dta$y2ana));
	#Firstly, to remove outlier in mid section
	y1<-predict(gam(y2ana~s(ts,bs='cc'),data=dta[ii_av,]),newdata=dta);
	yr<-as.vector((y1-dta$y2ana))
	q<-getOutliers(yr)
	iqq<-c(q$iRight,q$iLeft)
	ii_na_tol<-c(ii_na,iqq)
	dta$y2ana[ii_na_tol]<-NA
	ii_na<-which(is.na(dta$y2ana));
	ii_av<-which(!is.na(dta$y2ana));	
	dta$y2ana[ii_na]<-predict(gam(y2ana~s(ts,bs='tp'),data=dta[ii_av,]),newdata=dta[ii_na,]);
	aa<-quantile(y2ana,c(0.05,0.95),na.rm=T);
	yd<-which((y2ana-dta$y2ana)>tol*(aa[2]-aa[1]));
	dta$y2ana[yd]<-y2ana[yd];
	#Secondly, to remove outlier in side parts
	y1<-predict(gam(y2ana~s(ts,bs='cp'),data=dta[ii_av,]),newdata=dta);	
	yr<-as.vector((y1-dta$y2ana)/dta$y2ana)
	q<-getOutliers(yr)
	iqq<-c(q$iRight,q$iLeft)
	ii_na_tol<-c(iqq)
	dta$y2ana[ii_na_tol]<-NA
	ii_na<-which(is.na(dta$y2ana));
	ii_av<-which(!is.na(dta$y2ana));	
	dta$y2ana[ii_na]<-predict(gam(y2ana~s(ts,bs='tp'),data=dta[ii_av,]),newdata=dta[ii_na,]);
	dta$y2ana[yd]<-y2ana[yd];

	# plot(y2ana,pch=12,col=3,ylim=c(0,10))
    # points(dta$y2ana,pch=19)	
	dta$y2ana
}


rv_outliers_w<-function(y2ana,ts,tswn=c(lft,rgt),tol=0.05){
	dta<-data.frame(y2ana,ts);
	ii_av<-which(!is.na(dta$y2ana));
	ii_na<-which(is.na(dta$y2ana));
	#Firstly, to remove outlier in mid section
	y1<-predict(gam(y2ana~s(ts,bs='cc'),data=dta[ii_av,]),newdata=dta);
	yr<-as.vector((y1-dta$y2ana))
	q<-getOutliers(yr)
	iqq<-c(q$iRight,q$iLeft)
	ii_na_tol<-c(ii_na,iqq)
	dta$y2ana[ii_na_tol]<-NA
	ii_na<-which(is.na(dta$y2ana));
	ii_av<-which(!is.na(dta$y2ana));	
	dta$y2ana[ii_na]<-predict(gam(y2ana~s(ts,bs='tp'),data=dta[ii_av,]),newdata=dta[ii_na,]);
	aa<-quantile(y2ana,c(0.05,0.95),na.rm=T);
	# yd<-which((y2ana-dta$y2ana)>tol*(aa[2]-aa[1]));
	# dta$y2ana[yd]<-y2ana[yd];
	#Secondly, to remove outlier in side parts
	y1<-predict(gam(y2ana~s(ts,bs='cp'),data=dta[ii_av,]),newdata=dta);	
	yr<-as.vector((y1-dta$y2ana)/dta$y2ana)
	q<-getOutliers(yr)
	iqq<-c(q$iRight,q$iLeft)
	ii_na_tol<-c(iqq)
	dta$y2ana[ii_na_tol]<-NA
	ii_na<-which(is.na(dta$y2ana));
	ii_av<-which(!is.na(dta$y2ana));	
	dta$y2ana[ii_na]<-predict(gam(y2ana~s(ts,bs='tp'),data=dta[ii_av,]),newdata=dta[ii_na,]);
	yd<-which(((y2ana-dta$y2ana)>tol*(aa[2]-aa[1]))&(ts>=tswn[1])&(ts<=tswn[2]));
	dta$y2ana[yd]<-y2ana[yd];
	# plot(y2ana,pch=12,col=3,ylim=c(0,10))
    # points(dta$y2ana,pch=19)	
	dta$y2ana
}


pickup_data<-function(vv){
	options("stringsAsFactors"=FALSE);
	latlon<-as.data.frame(t(matrix(as.numeric(unlist(strsplit(gsub('Samp1Line1','',gsub("Lat",'',vv[,9])),'Lon'))),nrow=2)));
	names(latlon)<-c('lat','lon');
	
	cid<-as.data.frame(t(matrix(as.character(unlist(strsplit(as.character(vv[,6]),".",fixed=T))),nrow=6)),stringsAsFactors = FALSE);
	names(cid)<-c("Product","AYearDoy","tilt","version","AcDate","Bands");
	
	elev<-vv[,5];names(elev)<-'elev';
	
	val0<-vv[,11];names(val0)<-'val0';
	
	year=as.numeric(substr(vv[,8],2,5))
	doy<-as.numeric(substr(vv[,8],6,8))
	
	ydID<-cid[,1];ydID[]<-'XXXXXXXXXXX';names(ydID)<-'ydID'
	fdat<-data.frame(ydID,latlon,elev,year,doy,cid,val0,stringsAsFactors=FALSE);
	fdat
}

sumfun <- function(x, ...) {
    c(m = mean(x, ...), md = median(x, ...), v = sd(x, ...), l = length(x), q05 = quantile(x, 0.05, ...), q95 = quantile(x, 0.95, ...))
}

split.col<-function(fdat){
	split<-fdat$year+fdat$doy/400;
	fdat<-data.frame(fdat,split);
	bands<-levels(as.factor(fdat$Bands));
	v1<-fdat[(fdat$Bands==bands[1]),];
	names(v1)[13]<-bands[1]
	for (i in 2:length(bands)){
		v2<-fdat[(fdat$Bands==bands[i]),13:14];
		names(v2)[1]<-bands[i];
		v1=merge(v1,v2,by.x='split',by.y='split')
	}
	v1[,'split']<-NULL;v1[,'Bands']<-NULL;
	v1
}

options("stringsAsFactors"=FALSE) 
file_site<-'cprj_site.csv'
file2save<-'MOD09A1_Data.csv'
flst<-dir(pattern='.asc',path='.//MOD09A1_IndividualPixel',full.names = T)

vv<-read.csv(flst[1],header=FALSE)
vv<-pickup_data(vv)
fnx<-as.character(unlist(strsplit(toupper(basename(flst[1])),split='END'))[1])
latlon<-as.data.frame(t(matrix(as.numeric(unlist(strsplit(gsub("LAT",'',fnx),split='LON'))),nrow=2)));
names(latlon)<-c('lat','lon');

fdat<-split.col(vv)
fdat$lon<-latlon$lon;fdat$lat<-latlon$lat;
write.table(fdat,file2save,sep=',',row.names=F)
for (i in 2:length(flst)){
	vv<-read.csv(flst[i],header=FALSE)
	vv<-pickup_data(vv)
	fdat<-split.col(vv)
	
	fnx<-as.character(unlist(strsplit(toupper(basename(flst[i])),split='END'))[1])
#	print(i);print(fnx)
#	lat<-as.numeric(substr(fnx,4,10));lon<-as.numeric(substr(fnx,15,22))
	latlon<-as.data.frame(t(matrix(as.numeric(unlist(strsplit(gsub("LAT",'',fnx),split='LON'))),nrow=2)));
	names(latlon)<-c('lat','lon');
#	print(i);print(latlon)
    assign('fdat$lon',latlon$lon);assign('fdat$lat',latlon$lat);
 #   print(c(fdat$lon[1],fdat$lat[1]))
	write.table(fdat,file2save,sep=',',append=T,row.names=FALSE,col.names=FALSE)
	remove(latlon)
}
#ydID,ydID2,Longitude,Latitude,Sample_end_latest,SS
site<-read.csv(file_site)
names(site)[names(site) == "Longitude"] <- "lon"
names(site)[names(site) == "Latitude"] <- "lat"
names(site)[names(site) == "Sample_end_latest"] <- "end.date"

dat<-read.csv(file2save)
l.f<-nrow(site)
for (i in 1:l.f){
	ii<-which((dat$lon %in% site$lon[i])&(dat$lat %in% site$lat))
	if(length(ii)>0){dat$ydID[ii]<-site$ydID[i];print(c(i,dat$lon[ii[1]],dat$lat[ii[1]],site$ydID[i],site$lon[i],site$lat[i]))}
	rm(ii)
}
write.table(dat,file2save,sep=',',row.names=F)

