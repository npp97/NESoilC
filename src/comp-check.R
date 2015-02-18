# TODO: Add comment
# 
# Author: lenovo
###############################################################################
#keys<-c("KeyCode=","DateBegin=","DateEnd=","Project=","Location=","KeyNo=","SpeciesName=","Unit=")
#
#inp <- readLines(fname, ok = TRUE, warn = FALSE)
#header.begin <- grep("^HEADER:$", inp)
#header.end <- grep("^DATA:(Tree|Single)$", inp)
#n <- length(header.end)

require(dplR)
require(doBy)
require(stringr)
dict_sp<-read.csv('dict_sp.csv',header=F,colClasses=c('character'))
dict_site<-read.csv('dict_site.csv',header=F,colClasses=c('character'))
flst.fh<-dir(pattern='.fh',recursive=T,full.names = TRUE)
flst.rwl<-dir(pattern='.rwl',recursive=T,full.names = TRUE)
flst.all<-c(flst.fh,flst.rwl)

flst1<-basename(flst.all)
f.ext<-substr(flst1,nchar(flst1)-2,nchar(flst1))

flst<-sub("^([^.]*).*", "\\1", flst1)

site<-matrix(NA,ncol=1,nrow=length(flst))
mtt<-ydIDt<-mtd<-ydIDd<-age<-sp<-site
diam<-height<-mtt;

# tree species
for (i in 1:nrow(dict_sp)){ 
	ii<-grep(dict_sp$V1[i],substr(flst,5,nchar(flst)))
	if (length(ii)>0){sp[ii]<-dict_sp$V1[i];ii<-NULL;}
}

#mountain name 
for (j in 1:nrow(dict_site)){
	ii1<-grep(dict_site$V1[j],substr(flst,1,5))
	if (length(ii1)>0) {site[ii1]<-dict_site$V1[j];ii1<-NULL;}
}

#ydID&age
for (i in 1:length(flst)){
	print(i);print(flst.all[i])
	a<-read.rwl(flst.all[i])
	age[i]<-nrow(a)[1]
	year_sampled<-row.names(a)[nrow(a)]
	if(toupper(site[i])=='CBS'){b<-unlist(strsplit(flst1[i],split='-'))[1];b<-str_trim(b);ydIDd[i]<-paste('CBSSYF',year_sampled,substr(b,nchar(b)-2,nchar(b)),sep='')}
	if(toupper(site[i])=='WDS'){b<-unlist(strsplit(flst1[i],split='-'))[1];b<-str_trim(b);ydIDd[i]<-paste('WDSSYF',year_sampled,substr(b,nchar(b)-2,nchar(b)),sep='')}
	if(toupper(site[i])=='LGS'){b<-unlist(strsplit(flst1[i],split='-'))[1];b<-str_trim(b);ydIDd[i]<-paste('LGSSYF',year_sampled,substr(b,nchar(b)-2,nchar(b)),sep='')}
	if(toupper(site[i])=='LX'){ydIDd[i]<-paste('LXSDYF',year_sampled,substr(flst1[i],3,5),sep='')}
	if(toupper(site[i])=='LYL'){ydIDd[i]<-paste('LYLSYF',year_sampled,substr(flst1[i],12,14),sep='')}
	if(toupper(site[i])=='LYLS'){ydIDd[i]<-paste('LYLSYF',year_sampled,substr(flst1[i],5,7),sep='')}
	if(toupper(site[i])=='ZGCL'){ydIDd[i]<-paste('ZGCLYF',year_sampled,substr(flst1[i],5,7),sep='')}
	if(toupper(site[i])=='HDL'){ydIDd[i]<-paste('HDLSYF',year_sampled,substr(flst1[i],4,6),sep='')}
	if(toupper(site[i])=='HDLS'){ydIDd[i]<-paste('HDLSYF',year_sampled,substr(flst1[i],5,7),sep='')}
	if(toupper(site[i])=='HDLSY'){ydIDd[i]<-paste('HDLSYF',year_sampled,substr(flst1[i],13,15),sep='')} #11-13 or 13-15
	if(toupper(site[i])=='XXAL'){b<-unlist(strsplit(flst1[i],split='-'))[1];b<-str_trim(b);ydIDd[i]<-paste('XXALYF',year_sampled,substr(b,nchar(b)-2,nchar(b)),sep='')}
	if(toupper(site[i])=='DXAL'){b<-unlist(strsplit(flst1[i],split='-'))[1];b<-str_trim(b);ydIDd[i]<-paste('DXALYF',year_sampled,substr(b,nchar(b)-2,nchar(b)),sep='')}
	if(toupper(site[i])=='FC'){ydIDd[i]<-paste('LGSSYF',year_sampled,substr(flst1[i],3,5),sep='')}
	if(toupper(site[i])=='LTDZ'){ydIDd[i]<-paste('LGSSYF',year_sampled,substr(flst1[i],6,8),sep='')}
	if(toupper(site[i])=='FS'){ydIDd[i]<-paste('LGSSYF',year_sampled,substr(flst1[i],3,5),sep='')}
	if(toupper(site[i])=='TL'){ydIDd[i]<-paste('LGSSYF',year_sampled,substr(flst1[i],3,5),sep='')}	
}

grep('LYL',site)->ii
site[ii]<-'LYL'
grep('HDL',site)->ii
site[ii]<-'HDL'
grep('FS',site)->ii
site[ii]<-'LGS'
grep('FC',site)->ii
site[ii]<-'LGS'
grep('TL',site)->ii
site[ii]<-'LGS'
grep('LTDZ',site)->ii
site[ii]<-'LGS'



#Extract Diam and Height from strings like filename and KeyCode/Project
ii<-grep(pattern="D\\d{1,}|D\\d{1,}.\\d{1,}",flst1) #find the filename containing Height and Diameter
regmatches(flst1[ii],regexpr("D\\d{1,}|D\\d{1,}.\\d{1,}",flst1[ii]))->D     #got the charater
sub("H\\d{1,}|H\\d{1,}.\\d{1,}","",D)->D
diam[ii]<-as.numeric(substr(D,2,nchar(D)))

ii<-grep(pattern="H\\d{1,}|H\\d{1,}.\\d{1,}",flst1) #find the filename containing Height and Diameter
regmatches(flst1[ii],regexpr("H\\d{1,}|H\\d{1,}.\\d{1,}",flst1[ii]))->H
height[ii]<-as.numeric(substr(H,2,nchar(H)))

for (i in 1:length(f.ext)){
	if(toupper(f.ext[i])=='.FH'){
		inp <- readLines(flst.all[i], ok = TRUE, warn = FALSE)
		ii<-grep(pattern="H\\d{1,}|H\\d{1,}.\\d{1,}",inp)
		if(length(ii)>0){
			H<-regmatches(inp[ii],regexpr("H\\d{1,}|H\\d{1,}.\\d{1,}",inp[ii]));
			height[i]<-as.numeric(substr(H,2,nchar(H)))
		}
		ii<-grep(pattern="D\\d{1,}|D\\d{1,}.\\d{1,}",inp)
		if(length(ii)>0){
			D<-regmatches(inp[ii],regexpr("D\\d{1,}|D\\d{1,}.\\d{1,}",inp[ii]))
			diam[i]<-as.numeric(substr(D,2,nchar(D)))
		}	
	}
}

dinfo<-data.frame(ydIDd,age,site,sp,diam,height)
dinfo.id.m<-summaryBy(age~ydIDd,data=dinfo,FUN=mean)
dinfo.id.mx<-summaryBy(age~ydIDd,data=dinfo,FUN=max)

dinfo.id.s<-data.frame(dinfo.id.m,max=dinfo.id.mx[,2])


#Extract some information from
site_info1920<-read.csv('plot_info_1920.csv')
s2<-matrix(NA,ncol=1,nrow=nrow(site_info1920))
for (j in 1:nrow(dict_site)){
	ii1<-grep(dict_site$V1[j],substr(site_info1920$ydID,1,5))
	if (length(ii1)>0) {s2[ii1]<-dict_site$V1[j];ii1<-NULL;}
}

grep('LYL',s2)->ii
s2[ii]<-'LYL'
grep('HDL',s2)->ii
s2[ii]<-'HDL'


s2<-paste(s2,substr(site_info1920$ydID,13,15),sep='')
site_info=data.frame(ydIDd=s2,site_info1920)

#combine info from denchronsequences and field investigation
all<-merge(site_info,dinfo,all=TRUE)  

write.csv(all,file='all_age.csv')

write.csv(dinfo,file='age.csv')

#for (i in 1:length(f.ext)){
#	if(toupper(f.ext[i])=='.FH'){
#		inp <- readLines(flst.all[i], ok = TRUE, warn = FALSE)
#		ii<-grep(pattern="KeyCode=",inp)
#		if(!(length(ii)>0)){
#		file.rename(flst.all[i],sub('.fh','.bak_fh',flst.all[i]))
#		con<-file(flst.all[i],'w')
#		writeLines(inp[1],con=con)
#		writeLines(paste('KeyCode=',flst1[i],sep=''),con=con)
#		writeLines(inp[2:length(inp)],con=con)
#		close(con)
#			print(flst.all[i]);
#		}
#		}	
#	}

###processing of XXAL and DXAL data
#dy<-read.csv('年轮-样地对应表.csv',colClasses=c(rep('character',18)))
#fn.full<-dir(pattern='.fh',recursive=T)
#fn.bn<-basename(fn.full)
#fn.dir<-dirname(fn.full)
#fn<-data.frame(fn.full,TRFN=fn.bn,fn.dir)
#
#dy$TRFN <- gsub("-", "", dy$TRFN)
#fn$TRFN <- gsub("-","",fn$TRFN)
#dy$TRFN <- gsub(" ", "\\1", dy$TRFN)
#fn$TRFN <- gsub(" ","\\1",fn$TRFN)
#fndy<-merge(dy,fn)
#
#new_name<-paste(as.character(fndy$fn.dir),as.character(fndy$filename.new),sep='/')
#file.rename(from=as.character(fndy$fn.full),to=new_name)
#

write.csv(cbind(flst,site,sp),file='site_sp.csv')

library(doBy)
mem<-summaryBy(diam+height+age~ydID,data=all,FUN=mean,na.rm=T)
mxx<-summaryBy(diam+height+age~ydID,data=all,FUN=max,na.rm=T)
mm<-merge(mem,mxx)

soc_a<-read.csv('D:\\东北\\data\\土壤\\SOC_A.csv') #A layer soil carbon density
all<-read.csv('D:\\东北\\data\\土壤\\S)
merge(soc_a,mm,)


