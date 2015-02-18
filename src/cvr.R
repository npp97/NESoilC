# TODO: Add comment
# 
# Author: lenovo
###############################################################################


grpdist<-function(X){ #X is the matrix
	require(cluster)
	gr<-as.fata.frame(as.factor(X))
	distgr<-daisy(gr,'gower')
	distgr
}

require(SoDA);
require(FD);
require(cluster);

data<-read.csv('mbabs.csv')

# Convert latitude and longitude form into the corresponding coordinates in X (east-west) and Y (north-south) distances,
# the origin for the desired coordinates is the minimum values for the supplied latitude and longitude coordinates.

env_bak<-data[,c('X','Y','elev','aspect_d','P_mm','Ta_C')];
row.names(env_bak)<-data$ydID2

w0=c(1,1,1,1,3,3);
ii<-which(data$SX=='XXA')
env<-env_bak[ii,1:6]
w0=w0[1:4]
#XY<-geoXY(env$lat,env$long)
#env$lat<-XY[,1];env$long<-XY[,2];
env<-as.matrix(env);

#using cmean

# Calculate the Gower index (1983)


g<-gowdis(env,w0)
g<-dist(env)
# cluster
env.ch.ward<-hclust(g, method='ward');

# search the cutting level using Silhouette width
Nn<-nrow(env)
kk<-floor((Nn-1)/10)
kkm<-20:2:kk
asw<-numeric(length(kkm))
for (kj in 1:length(kkm)){
	sil<-silhouette(cutree(env.ch.ward,k=kkm[kj]),g)
	asw[kj]<-summary(sil)$avg.width
}
ks.best<-kkm[which.max(asw)]
ct<-cutree(env.ch.ward,k=ks.best)
ect=data.frame(env,ct,max_spec=data$max_spec[ii],max_spec2=data$max2_spec[ii])

plot(ect$X,ect$Y)
text(ect$X,ect$Y,label=ct)

# Search the cutting level using Mantel-optimal method
kt<-data.frame(k=1:nrow(env),r=0)
for (i in 2:(Nn-1)){
	gr<-cutree(env.ch.ward,i)
	distgr<-grpdist(gr)
	mt<-cor(g,distgr,method='pearson')
	kt[i,2]<-mt
}

km.best<-which.max(kt$r)

# Get the final tree
k.best<-which.max(ks.best,km.best)
env.ct<-cutree(env.ch.ward,k=k.best)

