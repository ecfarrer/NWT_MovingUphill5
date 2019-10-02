# Testing whether taking a subset of low and high with the same range of pca1 will give similar results (i.e. low has more complexity than high)

# I will start with the clr data sets from lo me hi. I don't think it is necessary to redo the clr claculations with the subst datset, b/c clr is calculated within sample (not between samples). the only difference would be that some taxa would be removed if they are not present in the subset plots, so it would only affect the number of 0 taxa

hmscYlo3
hmscYme3
hmscYhi3

ind<-which(comm.bio$lomehi=="lo")
hmscXlo<-data.frame(snowdepth=comm.bio$snowdepth,pH=comm.bio$pH,moisture=comm.bio$moisture,cvsnow=comm.bio$cvsnow,pca1=comm.bio$pca1)[ind,]

ind<-which(comm.bio$lomehi=="me")
hmscXme<-data.frame(snowdepth=comm.bio$snowdepth,pH=comm.bio$pH,moisture=comm.bio$moisture,cvsnow=comm.bio$cvsnow,pca1=comm.bio$pca1)[ind,]

ind<-which(comm.bio$lomehi=="hi")
hmscXhi<-data.frame(snowdepth=comm.bio$snowdepth,pH=comm.bio$pH,moisture=comm.bio$moisture,cvsnow=comm.bio$cvsnow,pca1=comm.bio$pca1)[ind,]

range(hmscXlo$pca1)
range(hmscXme$pca1)
range(hmscXhi$pca1)

#need a 0.384 pca1 window from high
temp<-sort(hmscXhi$pca1)
temp2<-rep(NA,10)
for (i in 1:13){
  a<-temp[i]+0.384
  temp2[i]<-length(which(temp<=a))-(i-1)
}
temp2
#take pca1 from 0.01448386 to 0.3984839, or 0.1074467 to 0.4914467

#lo
ind<-which(comm.bio$lomehi=="lo")
hmscXlo<-data.frame(snowdepth=comm.bio$snowdepth,pH=comm.bio$pH,moisture=comm.bio$moisture,cvsnow=comm.bio$cvsnow,pca1=comm.bio$pca1)[ind,]
rownames(hmscXlo)<-comm.bio$X.SampleID[ind]

which(hmscXlo$pca1<(-0.7386601)|(hmscXlo$pca1<(-0.354780)&hmscXlo$pca1>(-0.354781)))
set.seed(6);ind<-c(sample(1:25,size=10),11,21) #make sure 11 and 21 are not included
hmscXlo2<-data.frame(snowdepth=hmscXlo$snowdepth,pH=hmscXlo$pH,moisture=hmscXlo$moisture,cvsnow=hmscXlo$cvsnow)[ind,]
rownames(hmscXlo2)<-rownames(hmscXlo)[ind]

#subset y 
ind<-which(rownames(hmscYlo3)%in%rownames(hmscXlo2))
hmscYlo4<-hmscYlo3[ind,]

#sort them the same
hmscYlo5<-hmscYlo4[order(rownames(hmscYlo4)),]
hmscXlo3<-hmscXlo2[order(rownames(hmscXlo2)),]




ind<-which(comm.bio$lomehi=="hi")
hmscXhi<-data.frame(snowdepth=comm.bio$snowdepth,pH=comm.bio$pH,moisture=comm.bio$moisture,cvsnow=comm.bio$cvsnow,pca1=comm.bio$pca1)[ind,]
hmiscDISThi<-data.frame(X=comm.bio$X,Y=comm.bio$Y,elevation=comm.bio$elevation)[ind,]
rownames(hmscXhi)<-comm.bio$X.SampleID[ind]
rownames(hmiscDISThi)<-comm.bio$X.SampleID[ind]
