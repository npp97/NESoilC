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
xxalp<-read.csv('alltrees.csv',na.strings='NULL',colClasses=c(rep('factor',3),rep('numeric',3),rep('character',7)))
# ydID                zwmcbh xyfID  xj  jj  sg  zxg    whq shl comm                                 recorder                   recorderDate modifyDate
#    1 CBSSYF080628001   紫椴    13     1.3 1.6 2.4 <NA>  营养期       中       <NA>                                 曹伟、黄祥童、于兴华、李岩、高燕、刘巍、        80628      00:00.0
xxalp$cab<-xxalp$xj*xxalp$xj*pi/4.0;
xxalp$xxjj<-xxalp$xj/xxalp$jj;

tree.spec<-as.character(unlist(levels(as.factor(xxalp$zwmcbh))));
site.lst<-as.character(unlist(levels(as.factor(xxalp$ydID))));

nsite<-length(site.lst);
nspec<-length(tree.spec);

#QA XJ-JJ
#
ooxxjj<-(xxalp$xxjj>1)
ooxxjj2<-(xxalp$xxjj<=1)

for (i in 1:nsite){
	for (j in 1:nspec){
		ooii<-(xxalp$ydID %in% site.lst[i])&(xxalp$zwmcbh %in% tree.spec[j])
		ii<-which(ooii&ooxxjj)		
		if ((length(ii)>0)) {
			jj<-which(ooii&ooxxjj2)
			if(length(jj)>0){
				xxalp$jj[ii]<-predict(xxalp$xj[jj],lm(xxalp$jj[jj]~xxalp$xj[jj]));}}
#			xxalp$xxjj[ii]<-median(xxalp$xxjj[jj]); xxalp$jj[ii]=xxalp$xj[ii]/median(xxalp$xxjj[jj])}	
	}
}

ctbl<-as.data.frame(matrix(ncol=2,nrow=nsite*nspec))
names(ctbl)<-c('ydID','zwmcbh')
cc<-merge(tree.spec,site.lst)
ctbl[,1]<-cc[,2]; ctbl[,2]<-cc[,1]; rm(cc);

#ctbl[,1:2]<-as.character(ctbl[,1:2]);

#To generate site X species table for important value calculation
table.topr<-as.data.frame(matrix(0,ncol=nspec,nrow=nsite));
names(table.topr)<-tree.spec;
row.names(table.topr)<-site.lst;

#tables: importance value; relative density; relative coverage; relative frequency.
tbl.impV<-tbl.rden<-tbl.rcvg<-tbl.rfrq<-table.topr; 

rden.1<-summaryBy(xj~ydID+zwmcbh,data=xxalp,FUN=length); #density  by plot and speces
rfrq.1<-summaryBy(xj~ydID+zwmcbh+xyfID,data=xxalp,FUN=length);
rfrq2.1<-summaryBy(xyfID~ydID+zwmcbh,data=rfrq.1,FUN=length); # frequency by plot and speces
rcvg.1<-summaryBy(cab~ydID+zwmcbh,data=xxalp,FUN=sum); # total cutting area by plot and speces

# --to reshape these calculation into real table------

#----relative density :: rden.2
rden.2<-merge(ctbl,rden.1,all=T);
names(rden.2)[3]<-'rden'

chick_m <- melt(rden.2, id=1:2, na.rm=FALSE)
rden.2<-acast(chick_m, ydID ~ zwmcbh, sum, na.rm=T) 

#----relative frequency :: rfrq
rfrq.2<-merge(ctbl,rfrq2.1,all=T);
names(rfrq.2)[3]<-'rfrq'
#-------
chick_m <- melt(rfrq.2, id=1:2, na.rm=FALSE)
rfrq.2<-acast(chick_m, ydID ~ zwmcbh, sum, na.rm=T) 
#-------relative coverage :: rcvg.2
rcvg.2<-merge(ctbl,rcvg.1,all=T);
names(rcvg.2)[3]<-'rcvg'

chick_m <- melt(rcvg.2, id=1:2, na.rm=FALSE)
rcvg.2<-acast(chick_m, ydID ~ zwmcbh, sum, na.rm=T) 

#---------normalization to calculate important values：rimp
rden.3<-rden.2/apply(rden.2,1,sum,na.rm=T);
rfrq.3<-rfrq.2/36;
rcvg.3<-rcvg.2/apply(rcvg.2,1,sum,na.rm=T);

rimp<-(rden.3+rfrq.3+rcvg.3)/3;

#---------Simpson Index:: simp_idx
simp_idx= 1- apply(rden.3*rden.3,1,sum);

#dominant species :: dom_spec
dom_spec<-as.data.frame(matrix(ncol=4,nrow=nrow(rimp)));
row.names(dom_spec)<-row.names(rimp);
names(dom_spec)<-c('max_spec','max','max2_spec','max2')
f2max=apply(rimp,1,fst2max);
f2maxi<-apply(rimp,1,fst2maxi);

dom_spec[,c(1,3)]<-colnames(rimp)[t(f2maxi)];
dom_spec[,c(2,4)]<-t(f2max);

#dominant species and plot number of XXAL 
#   白桦 			斑叶稠李                稠李   	臭冷杉     	春榆     	         风桦    		蒿柳     	黑桦		 红皮云杉    		 红松  		 胡桃楸    	 櫰槐     	糠椴 
#      65        1        1       15       16        3        1        6        9       28        3        1        1 
# 		裂叶榆		  落叶松  	 毛赤杨   	蒙古栎   		青楷槭   	色木槭    	 山杨   		水曲柳   		  榆树  		 樟子松   	  紫椴 
#       1       80        6       44        1        5       11       10        1       21        1 

win.graph(width=6,height=4)
par(mfrow=c(1,2))
plot(density(dom_spec[,2]),main='第一建群种优势度概率密度');grid();abline(v=mean(dom_spec[,2]))
plot(density(dom_spec[,4]),main='第二建群种优势度概率密度');grid();abline(v=mean(dom_spec[,4]))

#---------Then calculate the correlationship between xj and jj
