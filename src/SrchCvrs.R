# TODO: Add comment
# 
# Author: lenovo
###############################################################################
# Grpdist
grpdist<-function(X){ #X is the matrix
	require(cluster)
	gr<-as.fata.frame(as.factor(X))
	distgr<-daisy(gr,'gower')
	distgr
}

require(SoDA);
require(FD);
require(cluster);

data<-read.csv('all_data.csv')

# Convert latitude and longitude form into the corresponding coordinates in X (east-west) and Y (north-south) distances,
# the origin for the desired coordinates is the minimum values for the supplied latitude and longitude coordinates.

env_bak<-data[,c('lat','long','elev','slope','aspect','posi','soil_depth')];
row.names(env_bak)<-data$ydID

w0=c(1,1,1,1,3,3,5);
env<-env_bak
XY<-geoXY(env$lat,env$long)
env$lat<-XY[,1];env$long<-XY[,2];
env<-as.matrix(env);

# Calculate the Gower index (1983)

g<-gowdis(env,w)

# cluster
env.ch.ward<-hclust(g, method='ward');

# search the cutting level using Silhouette width
Nn<-nrow(env)
asw<-numeric(Nn)

for (k in 2:(Nn-1)){
	sil<-silhouette(cutree(env.ch.ward,k=k),g)
	asw[k]<-summary(sil)$avg.width
}

ks.best<-which.max(asw)

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


