# TODO: Add comment
# 
# Author: lenovo
###############################################################################
require(dplR)
require(pander)
require(doBy)
#----------------------
sumfun <- function(x, ...){
	c(m=mean(x, ...), md=median(x,...),v=sd(x, ...), l=length(x,...),mi=min(x,...),mx=max(x,...))
}
#--------------
loc.spe<-function(tmp){
	spe<-locc<-ag<-NA
	spp<-agrep('SpeciesName',tmp);
	locp<-agrep('Location',tmp);
	agp<-agrep('Length',tmp)
	if(spp) spe<-trim(unlist(strsplit(tmp[spp],'='))[2]);
	if(locp) locc<-trim(unlist(strsplit(unlist(strsplit(tmp[locp],'='))[2],'-'))[1]);
	if(agp) ag<-as.numeric(unlist(strsplit(tmp[agp],'='))[2])
	return(c(locc,spe,ag));
}
#--------------------------------------------------------------------------------------------------------------------------------------
#please set the data fold as the working fold. 

flst<-dir(pattern='*.fh', recursive = TRUE,ignore.case = TRUE)
	lflst<-length(flst)
	las<-as.data.frame(matrix(NA,ncol=5,nrow=lflst))
	names(las)<-c('TRFN','location','species','age','filename')	
	chr.rwl<-read.rwl(flst[1],format='auto')
	con <- file(flst[1], "r", blocking = FALSE)
		tmp<-readLines(con,10)
	close(con)
	las[1,5]<-flst[1]
	las[1,2:4]<-loc.spe(tmp);
	

for (i in 2:lflst){
	try({
	chr1.rwl <- read.rwl(flst[i],format='auto');
	mxm<-max(as.numeric(row.names(chr1.rwl)))
	mim<-min(as.numeric(row.names(chr1.rwl)))
	if((mxm<2012)&(mim>1730)){
		chr.rwl = combine.rwl(chr.rwl, chr1.rwl);
		con <- file(flst[i], "r", blocking = FALSE)
			tmp<-readLines(con,10)
		close(con)
		las[i,5]<-flst[i]
		las[i,2:4]<-loc.spe(tmp);
	}})
}

save(chr.rwl,file='chr_rwl.RData')
plot(chron(chr.rwl))

las$age<-as.numeric(las$age)
las$species<-as.factor(las$species)

require(doBy)
#summary(age~species,data=las)
summaryBy(age~species,data=las,FUN=sumfun)

diam_age<-chr.rwl
diam_age[]<-NA
coef_6reg<-as.data.frame(matrix(nrow=ncol(chr.rwl),ncol=8));
names(coef_6reg)<-c('file_name','diam','I(diam^2)','I(diam^3)','I(diam^4)','I(diam^5)','I(dima^6)','adj-R2')

for (i in 1:ncol(chr.rwl)){
	i.nna<-which(!is.na(chr.rwl[,i]))
	diam_age[1:length(i.nna),i]<-chr.rwl[i.nna,i]
	x<-chr.rwl[i.nna,i];
	y<-1:length(i.nna);
	m<-lm(y~x+I(x^2)+I(x^3)+I(x^4)+I(x^5)+I(x^6)-1);
	coef_6reg[i,1]<-flst[i];
	coef_6reg[i,2:7]<-coef(m);
	coef_6reg[i,8]<-summary(m)[9];
}

save(las,chr.rwl,coef_6reg,diam_age,file='Fst_save.RData')

site_tring<-read.csv('年轮-样地对应表.csv');
site_tring$ydID<-as.character(site_tring$ydID);
site_tring$TRFN<-as.character(site_tring$TRFN);

for (i in 1:nrow(las)){
	m<-unlist(strsplit(las$filename[i],'/'))
	las$TRFN[i]<-m[length(m)]
}

las$TRFN<-gsub('-','',las$TRFN)
site_tring$TRFN<-gsub('-','',site_tring$TRFN)
site_info<-merge(las,site_tring)


