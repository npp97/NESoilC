# TODO: Add comment
# 
# Author: lenovo
###############################################################################
require(doBy)
require(lattice)
sumfun <- function(x, ...){
	c(m=mean(x, ...), md=median(x,...),v=sd(x, ...), l=length(x))
}

setwd("C:\\Users\\lenovo\\workspace\\NESoilC\\data")

w_xxal<-read.csv('XXAL.csv');
head(w_xxal)
xxal<-w_xxal[,c(2:4,6,8,10,12,13:17)];
xxal$id<-1:nrow(xxal)
xxal$soil_Ccb<-(xxal$soil_C^0.25-1)/0.25;
head(xxal);
# >        lat     long elv A_thick AAB_tck soil_depth    Diam_C soil_C     abv_C
# 1 47.36000 128.1930 305      15      45         65 13.920833     11 2.8025015
#    BAR                   comm   dominant
# 1 0.25       兴安落叶松人工林 兴安落叶松

w_xxal<-read.csv('xxal_ful.csv');
head(w_xxal)

#3----------------------
xxalf<-w_xxal[,c(2:4,10,12,14,15)];
xxalf$id<-1:nrow(xxal)
xxalf$soil_Ccb<-(xxal$soil_C^0.25-1)/0.25;
head(xxalf);

rm(w_xxal);
save(xxal,xxalf,file='xxal.RData')

#-----------------------------XXAL_A


yy.sp<-as.numeric(cut(xxal$Diam_C,breaks=quantile(xxal$Diam_C,probs=seq(0.,1,.1),na.rm=T)))
summaryBy(Diam_C+abv_C+soil_C~yy.sp,data=data.frame(xxal,yy.sp),na.rm=T,FUN=sumfun)->sc.sf
sc.sf<-sc.sf[1:10,]
sc.sf$soil_C.dp<-sc.sf$soil_C.v/sqrt(sc.sf$soil_C.l+0.001);

require(Hmisc)

win.graph(height=6,width=7)
par(mfrow=c(2,2),mar=c(4,4,2,1))

plot(1,xlim=c(5,25),ylim=c(4,11),type='n',xlab='平均胸径(cm)',ylab='表层土壤碳密度')
errbar(x=sc.sf$Diam_C.m,y=sc.sf$soil_C.m,yplus=sc.sf$soil_C.m+sc.sf$soil_C.dp,yminus=sc.sf$soil_C.m-sc.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
axis(1);axis(2);grid()
bb<-(lm((soil_C.m)~Diam_C.m+I(Diam_C.m^2)+I(Diam_C.m^3),data=sc.sf))  # 2nd order regression Y = b0 + (b1 * t) + (b2 * t**2)
lines(sc.sf$Diam_C.m,(predict(bb)),lwd=2,col=2)
title(main='平均值')

plot(1,xlim=c(5,25),ylim=c(4,11),type='n',xlab='平均胸径(cm)',ylab='表层土壤碳密度')
errbar(x=sc.sf$Diam_C.md,y=sc.sf$soil_C.md,yplus=sc.sf$soil_C.md+sc.sf$soil_C.dp,yminus=sc.sf$soil_C.md-sc.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
axis(1);axis(2);grid()
bb<-(lm((soil_C.md)~Diam_C.md+I(Diam_C.md^2)+I(Diam_C.md^3),data=sc.sf))  # 2nd order regression Y = b0 + (b1 * t) + (b2 * t**2)
lines(sc.sf$Diam_C.md,(predict(bb)),lwd=2,col=2)
title(main='中值')

plot(1,xlim=c(1,7),ylim=c(4,11),type='n',xlab='地上生物量',ylab='表层土壤碳密度')
errbar(x=sc.sf$abv_C.m,y=sc.sf$soil_C.m,yplus=sc.sf$soil_C.m+sc.sf$soil_C.dp,yminus=sc.sf$soil_C.m-sc.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
axis(1);axis(2);grid()
bb<-(lm((soil_C.m)~abv_C.m+I(abv_C.m^2)+I(abv_C.m^3),data=sc.sf))  # 2nd order regression Y = b0 + (b1 * t) + (b2 * t**2)
lines(sc.sf$abv_C.m,(predict(bb)),lwd=2,col=3)


plot(1,xlim=c(1,7),ylim=c(4,11),type='n',xlab='地上生物量',ylab='表层土壤碳密度')
errbar(x=sc.sf$abv_C.md,y=sc.sf$soil_C.md,yplus=sc.sf$soil_C.md+sc.sf$soil_C.dp,yminus=sc.sf$soil_C.md-sc.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
axis(1);axis(2);grid()
bb<-(lm((soil_C.md)~abv_C.md+I(abv_C.md^2)+I(abv_C.md^3),data=sc.sf))  # 2nd order regression Y = b0 + (b1 * t) + (b2 * t**2)
lines(sc.sf$abv_C.md,(predict(bb)),lwd=2,col=3)


#bb1<-(lm(log(soil_C.m)~I(1/Diam_C.m),data=sc.sf)) #S curve:Y = e**(b0 + (b1/t)) or ln(Y) = b0 + (b1/t)
#lines(sc.sf$Diam_C.m,exp(predict(bb1)),lwd=2,col=3)
#
#yy.sp<-as.numeric(cut(xxal$abv_C,breaks=quantile(xxal$abvC,probs=seq(0.,1,.1),na.rm=T)))
#summaryBy(Diam_C+abv_C+soil_C~yy.sp,data=data.frame(xxal,yy.sp),na.rm=T,FUN=sumfun)->sc.sf
#sc.sf<-sc.sf[1:10,]
#sc.sf$soil_C.dp<-sc.sf$soil_C.v/sqrt(sc.sf$soil_C.l+0.001);
#plot(1,xlim=c(1,7),ylim=c(3,12),type='n',xlab='乔木碳储量',ylab='表层土壤碳密度')
#errbar(x=sc.sf$abv_C.m,y=sc.sf$soil_C.md,yplus=sc.sf$soil_C.md+sc.sf$soil_C.dp,yminus=sc.sf$soil_C.md-sc.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
#axis(1);axis(2);grid()
#
#
#
#yy1.sp<-as.numeric(cut(xxalf$abv_C,breaks=quantile(xxalf$abv_C,probs=seq(0.,1,.1),na.rm=T)))
#summaryBy(Diam_C+abv_C+soil_C~yy1.sp,data=data.frame(xxalf,yy1.sp),na.rm=T,FUN=sumfun)->sc.sf1
#sc.sf1<-sc.sf1[1:10,]
#sc.sf1$soil_C.dp<-sc.sf1$soil_C.v/sqrt(sc.sf1$soil_C.l+0.001);
#
#plot(1,xlim=c(1,7),ylim=c(3,12),type='n',xlab='乔木碳储量',ylab='表层土壤碳密度')
#errbar(x=sc.sf1$abv_C.m,y=sc.sf1$soil_C.m,yplus=sc.sf1$soil_C.m+sc.sf1$soil_C.dp,yminus=sc.sf1$soil_C.m-sc.sf1$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
#axis(1);axis(2);grid()

#---------------------------------------------XXAL_FULL---------

#xyplot(soil_C~Diam_C|dominant,data=xxalf,type = c("smooth",'p'))
#xyplot(soil_C~Diam_C,data=xxalf,type = c("smooth",'p'))
#
#xyplot(soil_C~abv_C|dominant,data=xxalf,type = c("smooth",'p'))
#xyplot(soil_C~abv_C,data=xxalf,type = c("smooth",'p'))
#

win.graph(height=6,width=7)
par(mfrow=c(2,2),mar=c(4,4,2,1))

yyf.sp<-as.numeric(cut(xxalf$Diam_C,breaks=quantile(xxalf$Diam_C,probs=seq(0.,1,.1),na.rm=T)))
summaryBy(Diam_C+abv_C+soil_C~yyf.sp,data=data.frame(xxalf,yyf.sp),na.rm=T,FUN=sumfun)->scf.sf
scf.sf<-scf.sf[1:10,]
scf.sf$soil_C.dp<-scf.sf$soil_C.v/sqrt(scf.sf$soil_C.l+0.001);
#
#require(Hmisc)
#win.graph(height=4,width=8)
#par(mfrow=c(1,2),mar=c(4,4,2,1))

plot(1,xlim=c(5,25),ylim=c(4,11),type='n',xlab='平均胸径(cm)',ylab='表层土壤碳密度')
errbar(x=scf.sf$Diam_C.m,y=scf.sf$soil_C.m,yplus=scf.sf$soil_C.m+scf.sf$soil_C.dp,yminus=scf.sf$soil_C.m-scf.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
axis(1);axis(2);grid()
bb<-(lm((soil_C.m)~Diam_C.m+I(Diam_C.m^2)+I(Diam_C.m^3),data=scf.sf))
lines(scf.sf$Diam_C.m,(predict(bb)),lwd=2,col=2)
title(main='平均值')

plot(1,xlim=c(5,25),ylim=c(4,11),type='n',xlab='平均胸径(cm)',ylab='表层土壤碳密度')
errbar(x=scf.sf$Diam_C.md,y=scf.sf$soil_C.md,yplus=scf.sf$soil_C.md+scf.sf$soil_C.dp,yminus=scf.sf$soil_C.md-scf.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
axis(1);axis(2);grid()
bb<-(lm((soil_C.md)~Diam_C.md+I(Diam_C.md^2)+I(Diam_C.md^3),data=scf.sf))
lines(scf.sf$Diam_C.md,(predict(bb)),lwd=2,col=2)
title(main='中间值')

plot(1,xlim=c(1,7),ylim=c(4,11),type='n',xlab='地上生物量',ylab='表层土壤碳密度')
errbar(x=scf.sf$abv_C.m,y=scf.sf$soil_C.m,yplus=scf.sf$soil_C.m+scf.sf$soil_C.dp,yminus=scf.sf$soil_C.m-scf.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
axis(1);axis(2);grid()
bb<-(lm((soil_C.m)~abv_C.m+I(abv_C.m^2)+I(abv_C.m^3),data=scf.sf))  # 2nd order regression Y = b0 + (b1 * t) + (b2 * t**2)
lines(scf.sf$abv_C.m,(predict(bb)),lwd=2,col=3)


plot(1,xlim=c(1,7),ylim=c(4,11),type='n',xlab='地上生物量',ylab='表层土壤碳密度')
errbar(x=scf.sf$abv_C.md,y=scf.sf$soil_C.md,yplus=scf.sf$soil_C.md+scf.sf$soil_C.dp,yminus=scf.sf$soil_C.md-scf.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
axis(1);axis(2);grid()
bb<-(lm((soil_C.md)~abv_C.md+I(abv_C.md^2)+I(abv_C.md^3),data=scf.sf))  # 2nd order regression Y = b0 + (b1 * t) + (b2 * t**2)
lines(scf.sf$abv_C.md,(predict(bb)),lwd=2,col=3)

#plot(1,xlim=c(1,7),ylim=c(4,11),type='n',xlab='乔木碳储量',ylab='表层土壤碳密度')
#errbar(x=sc.sf$abv_C.m,y=sc.sf$soil_C.md,yplus=sc.sf$soil_C.md+sc.sf$soil_C.dp,yminus=sc.sf$soil_C.md-sc.sf$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
#axis(1);axis(2);grid()
#
#
#yy1.sp<-as.numeric(cut(xxalf$abv_C,breaks=quantile(xxalf$abv_C,probs=seq(0.,1,.1),na.rm=T)))
#summaryBy(Diam_C+abv_C+soil_C~yy1.sp,data=data.frame(xxalf,yy1.sp),na.rm=T,FUN=sumfun)->sc.sf1
#sc.sf1<-sc.sf1[1:10,]
#sc.sf1$soil_C.dp<-sc.sf1$soil_C.v/sqrt(sc.sf1$soil_C.l+0.001);
#
#plot(1,xlim=c(1,7),ylim=c(3,12),type='n',xlab='乔木碳储量',ylab='表层土壤碳密度')
#errbar(x=sc.sf1$abv_C.m,y=sc.sf1$soil_C.m,yplus=sc.sf1$soil_C.m+sc.sf1$soil_C.dp,yminus=sc.sf1$soil_C.m-sc.sf1$soil_C.dp,cex=1.1,add=T,pch=19,cap=0);box()
#axis(1);axis(2);grid()
#
#


#bb1<-(lm(log(soil_C.m)~I(Diam_C.m),data=scf.sf)) #S curve:Y = e**(b0 + (b1/t)) or ln(Y) = b0 + (b1/t)
#lines(scf.sf$Diam_C.m,exp(predict(bb1)),lwd=2,col=3)
#

win.graph()
xyplot(soil_C~Diam_C|dominant,data=xxal,type = c("smooth",'p'))

win.graph()
xyplot(soil_C~Diam_C,data=xxal,type = c("smooth",'p'))

win.graph()
xyplot(soil_C~Diam_C|dominant,data=xxalf,type = c("smooth",'p'))

win.graph()
xyplot(soil_C~Diam_C,data=xxalf,type = c("smooth",'p'))

#xyplot(soil_C~abv_C|dominant,data=xxal,type = c("smooth",'p'))
#xyplot(soil_C~abv_C,data=xxal,type = c("smooth",'p'))
#


