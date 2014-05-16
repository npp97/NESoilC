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

loc.spe<-function(tmp){
	spe<-locc<-ag<-NA
	spp<-agrep('SpeciesName',tmp);
	locp<-agrep('Location',tmp);
	agp<-agrep('Length',tmp)
	if(spp) spe<-trim.spaces(unlist(strsplit(tmp[spp],'='))[2]);
	if(locp) locc<-trim.spaces(unlist(strsplit(unlist(strsplit(tmp[locp],'='))[2],'-'))[1]);
	if(agp) ag<-as.numeric(unlist(strsplit(tmp[agp],'='))[2])
	return(c(locc,spe,ag));
}

l.nna<-function(x){
	return(length(which(!is.na(x))));
}

to.chron.by.col<-function(x){
	nc<-names(x);
	lvl.nc<-levels(as.factor(nc));
	l.lnc<-length(lvl.nc);
	ch.p.site<-as.data.frame(matrix(NA,ncol=l.lnc,nrow=nrow(x)))
	names(ch.p.site)<-as.character(lvl.nc)
	for (i in 1:l.lnc){
		ii<-which(nc %in% lvl.nc[i])
		if (length(ii)<2) ch.p.site[,i]<-x[,ii]
		else ch.p.site[,i]<-chron(x[,ii])[1];
	}
	return(ch.p.site);
}

grw_crv<-function(pars,x,y){
	m<-pars[1]*x^pars[1]+pars[3];
	return(sum((m-y)^2))
}

#--------------------------------------------------------------------------------------------------------------------------------------
#please set the data fold as the working fold. 

#read raw tree-ring chronsequences
flst<-dir(pattern='*.fh', recursive = TRUE,ignore.case = TRUE)
	lflst<-length(flst)
#read site_info table	
	las<-as.data.frame(matrix(NA,ncol=5,nrow=lflst))
	names(las)<-c('TRFN','location','species','age','filename')	
#begin to read raw chronsequence and to build the cross-section of raw-chron and site-info
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
#		chr.rwl = combine.rwl(chr.rwl, chr1.rwl);
		con <- file(flst[i], "r", blocking = FALSE)
			tmp<-readLines(con,10)
		close(con)
		las[i,5]<-flst[i]
		las[i,2:4]<-loc.spe(tmp);
	}})
}

#save(chr.rwl,file='chr_rwl.RData')
#plot(chron(chr.rwl))
for (i in 1:nrow(las)){
	m<-unlist(strsplit(las$filename[i],'/'))
	las$TRFN[i]<-m[length(m)]
}

las$age<-as.numeric(las$age)
las$species<-as.factor(las$species)
las$TRFN<-gsub('-','',las$TRFN)

site_tring<-read.csv('年轮-样地对应表.csv');
site_tring$ydID<-as.character(site_tring$ydID);
site_tring$TRFN<-as.character(site_tring$TRFN);
site_tring$TRFN<-gsub('-','',site_tring$TRFN)
site_info<-merge(las,site_tring)

l.lst<-nrow(site_info)

chr.rwl<-read.rwl(site_info$filename[1],format='auto')
for (i in 2:l.lst){
	try({
				chr1.rwl <- read.rwl(site_info$filename[i],format='auto');
				mxm<-max(as.numeric(row.names(chr1.rwl)))
				mim<-min(as.numeric(row.names(chr1.rwl)))
				if((mxm<2012)&(mim>1730)){
					chr.rwl = combine.rwl(chr.rwl, chr1.rwl);
				}})
}

site_info$age0<-apply(chr.rwl,2,FUN=l.nna);
site_info$树种<-trim.spaces(site_info$树种);
site_info$ss<-paste(site_info$ydID,site_info$树种,sep='-');

#Real site info, grouped by ydID and 树种
site_info_by_ss<-site_info[unique(site_info$ss),];

#-----------------------------------------------
# site_info :: "TRFN","location","species","age","filename","ydID","树种","qlName","region","latitudeN", "longitudeE","altitude","pd","px","pw","zyyssz","ysgm","yscb","备注","age0"  
#-----------------------------------------------
names(chr.rwl)<-paste(site_info$ydID,site_info$树种,sep='-')
chron.site<-to.chron.by.col(chr.rwl)

require(doBy)
#summary(age~species,data=las)
summaryBy(age0~树种,data=site_info,FUN=sumfun)
summaryBy(age0~ydID,data=site_info,FUN=sumfun)
summaryBy(age0~ydID+树种,data=site_info,FUN=sumfun)
plot(density(summaryBy(age0~ydID,data=site_info,FUN=sumfun)$age0.mx),main='大小兴安岭调查群落最大年龄分布')

#to build the age-diam curves using A=aD^b
diam_age<-chr.rwl
diam_age[]<-NA
coef_6reg<-as.data.frame(matrix(nrow=ncol(chr.rwl),ncol=15));
names(coef_6reg)<-c('file_name','avg_grw_rate_05','avg_grw_rate_all','modTpe0','a0','b0','c0','modTpe1','a1','b1','c1','n0','n1','sigma0','sigma1')
start<-c(a=0.5,b=0.5,c=1)
start<-c(a=0.5,b=0.5)
resids<-chr.rwl;
resids0<-resids;

for (i in 1:ncol(chr.rwl)){
	i.nna<-which(!is.na(chr.rwl[,i]))
	diam_age[1:length(i.nna),i]<-chr.rwl[i.nna,i]

	coef_6reg[i,1]<-flst[i];
	coef_6reg[i,2]<-mean(chr.rwl[i.nna[min(5,length(i.nna))],i]);
	coef_6reg[i,3]<-mean(chr.rwl[i.nna,i]);
	coef_6reg[i,13]<-length(i.nna);	
	coef_6reg[i,12]<-min(5,length(i.nna));
	
	x<-cumsum(chr.rwl[i.nna,i])*2; #convert ring-width into diameter
	y<-1:length(i.nna);
	
	st<-try(nls(y~a*x^b,start=start))
	
	if(class(st)!='try-error'){
		m1<-nls(y~a*x^b,start=start);
		resids[i.nna,i]<-resid(m1);
		coef_6reg[i,9:10]<-coef(m1);
		coef_6reg[i,15]<-summary(m1)$sigma;
		coef_6reg[i,8]<-'expM';		
	}
	else{
		m1<-lm(y~x+I(x^2));
		resids[i.nna,i]<-resid(m1);
		coef_6reg[i,9:10]<-coef(m1);
		coef_6reg[i,15]<-summary(m1)$sigma;
		coef_6reg[i,8]<-'LM';
	}


	x<-x[1:min(5,length(i.nna))];
	y<-y[1:min(5,length(i.nna))];
	
	st<-try(nls(y~a*x^b,start=start))
	if(class(st)!='try-error'){
		m1<-nls(y~a*x^b,start=start);
		resids0[i.nna[1:min(5,length(i.nna))],i]<-resid(m1);
		coef_6reg[i,5:6]<-coef(m1);
		coef_6reg[i,14]<-summary(m1)$sigma;
		coef_6reg[i,4]<-'expM';		
	}
	else{
		m1<-lm(y~x+I(x^2));
		resids0[i.nna[1:min(5,length(i.nna))],i]<-resid(m1);
		coef_6reg[i,5:6]<-coef(m1);
		coef_6reg[i,14]<-summary(m1)$sigma;
		coef_6reg[i,4]<-'LM';
	}	
}

diam.age.chron.site<-to.chron.by.col(diam_age)
#Now, to calculate the site-specie Age-Diam curves (1) Age=a*Diam^b or (2) Age=a+b*Diam+c*Diam^2
coef.diam.reg.site<-as.data.frame(matrix(ncol=15,nrow=ncol(chron.site)));
names(coef.diam.reg.site)<-c('site_specie','avg_grw_rate_05','avg_grw_rate_all','modTpe0','a0','b0','c0','modTpe1','a1','b1','c1','n0','n1','sigma0','sigma1');
start<-c(a=0.5,b=0.5)
resids.ss<-chron.site;
resids0.ss<-resids.ss;

for (i in 1:ncol(chron.site)){
	i.nna<-which(!is.na(chron.site[,i]))
	
	coef.diam.reg.site[i,1]<-names(chron.site)[i];
	coef.diam.reg.site[i,2]<-mean(chron.site[i.nna[min(5,length(i.nna))],i]);
	coef.diam.reg.site[i,3]<-mean(chron.site[i.nna,i]);
	coef.diam.reg.site[i,13]<-length(i.nna);	
	coef.diam.reg.site[i,12]<-min(5,length(i.nna));
	
	x<-cumsum(chron.site[i.nna,i])*2; #convert ring-width into diameter
	y<-1:length(i.nna);
	
	st<-try(nls(y~a*x^b,start=start))
	
	if(class(st)!='try-error'){
		m1<-nls(y~a*x^b,start=start);
		resids[i.nna,i]<-resid(m1);
		coef.diam.reg.site[i,9:10]<-coef(m1);
		coef.diam.reg.site[i,15]<-summary(m1)$sigma;
		coef.diam.reg.site[i,8]<-'expM';		
	}
	else{
		m1<-lm(y~x+I(x^2));
		resids[i.nna,i]<-resid(m1);
		coef.diam.reg.site[i,9:10]<-coef(m1);
		coef.diam.reg.site[i,15]<-summary(m1)$sigma;
		coef.diam.reg.site[i,8]<-'LM';
	}

	
	x<-x[1:min(5,length(i.nna))];
	y<-y[1:min(5,length(i.nna))];
	
	st<-try(nls(y~a*x^b,start=start))
	if(class(st)!='try-error'){
		m1<-nls(y~a*x^b,start=start);
		resids0[i.nna[1:min(5,length(i.nna))],i]<-resid(m1);
		coef.diam.reg.site[i,5:6]<-coef(m1);
		coef.diam.reg.site[i,14]<-summary(m1)$sigma;
		coef.diam.reg.site[i,4]<-'expM';		
	}
	else{
		m1<-lm(y~x+I(x^2));
		resids0[i.nna[1:min(5,length(i.nna))],i]<-resid(m1);
		coef.diam.reg.site[i,5:6]<-coef(m1);
		coef.diam.reg.site[i,14]<-summary(m1)$sigma;
		coef.diam.reg.site[i,4]<-'LM';
	}	
}

save(chr.rwl,diam_age,              #two kind of rqw chronsequences
	chron.site,diam.age.chron.site, #chronsequences grouped by ydID-树种
	coef_6reg,coef.diam.reg.site,   #Age-Diam curves of each tree-ring series and mean chrons by ydID-树种 
	las,site_info,site_info_by_ss,  #Information of each site, grouped by chron and ydID-树种 respectatively.
	resids,resids0,file='Fst_save.RData') # residual values of fitted models

#todo 
#	to obtain relationship between xj and jj for dominant species of each site
#   objectives : (1) interception ; (2) slope 
#   input: site_info_ss :: to provide 
#   
Enter file contents herealltrees<-read.csv('alltrees.csv',na.strings='NULL',colClasses=c(rep('factor',3),rep('numeric',3),rep('character',7)))
# ydID                zwmcbh xyfID  xj  jj  sg  zxg    whq shl comm                                 recorder                   recorderDate modifyDate
#    1 CBSSYF080628001   紫椴    13     1.3 1.6 2.4 <NA>  营养期       中       <NA>                                 曹伟、黄祥童、于兴华、李岩、高燕、刘巍、        80628      00:00.0
with(alltrees,alltrees$xxjj<-xj/jj);