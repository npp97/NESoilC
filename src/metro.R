# TODO: to Generate annaual temperature and precipation using GAM, driven by lat,long and elev, and IDW model. the GAM model produce the smooth surface and IDW model generated residuals 
# 
# Author: Junhui@Nov,2014
###############################################################################


require(SoDA);
require(mgcv);
require(gstat);
require(sp)
require(geoR)

met0<-read.csv('met_ne.csv');
met_p<-read.csv('biomass.csv');

XY0<-geoXY(met0$lat_d,met0$long_d,lat0=38.9, lon0=119,unit=1000)
XY_p<-geoXY(met_p$lat,met_p$long,lat0=38.9, lon0=119,unit=1000)

met0<-data.frame(met0,XY0);
rs.grid<-expand.grid(X=seq(-160,1062,1),Y=seq(110,1604,1))

met_p<-data.frame(met_p,XY_p);
ii<-which((XY_p[,1]<50)&(XY_p[,2]<500))
met_p[-ii,]->met_p

#predict annual precipitation using lat, long and elev
p_model<-gam(p_mm~s(X)+s(Y)+s(elev),data=met0,select=T);
p_prd<-as.data.frame(predict.gam(p_model,met_p));
names(p_prd)<-"p_prd"
rs_p<-resid(p_model);

ta_model<-gam(ta_c~s(X)+s(Y)+s(elev),data=met0,select=T);
ta_prd<-as.data.frame(predict.gam(ta_model,met_p));
names(ta_prd)<-"ta_prd";

rs_ta<-resid(ta_model);

met0<-data.frame(met0,rs_p=rs_p,rs_ta=rs_ta)
#met_p<-data.frame(met_p,p_prd=p_prd,ta_prd=ta_prd)

coordinates(met0) = ~X+Y
gridded(rs.grid) = ~X+Y

pr.idw<-as.data.frame(idw(rs_p~1,met0,rs.grid))
ta.idw<-as.data.frame(idw(rs_ta~1,met0,rs.grid))

np<-nrow(met_p)
rs_ta<-rs_p<-matrix(NA,np)
rs.grid<-as.data.frame(rs.grid)
for (i in 1:np){ #predict points
	ii<-which((rs.grid$x<(1.2*met_p$X[i]))&(rs.grid$x>(0.8*met_p$X[i]))&(rs.grid$y<(1.2*met_p$Y[i]))&(rs.grid$y>(0.8*met_p$Y[i])));
	nr<-length(ii);
	if(nr>0){
	dis0<-matrix(NA,nr,2);	
	dis0[,1]<-met_p$X[i];
	dis0[,2]<-met_p$Y[i];
	
	dis1<-sqrt((dis0[,1]-rs.grid[ii,1])*(dis0[,1]-rs.grid[ii,1])+(dis0[,2]-rs.grid[ii,2])*(dis0[,2]-rs.grid[ii,2]));
	ij<-which.min(dis1);
	rs_p[i]<-pr.idw[ii[ij],3];
	rs_ta[i]<-ta.idw[ii[ij],3];
	print(c(i,dis1[ij],rs_p[i],rs_ta[i]))
	}
	else print(i)
 
}


p_fin<-p_prd+rs_p;
ta_fin<-ta_prd+rs_ta;
names(p_fin)<-"p_fin";
names(ta_fin)<-"ta_fin";

met_f<-data.frame(met_p,p_prd,ta_prd,p_fin,ta_fin)

write.csv(met_f,file='met_p.csv',row.names=F)

#meteorl



