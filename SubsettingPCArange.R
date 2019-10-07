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
#range is 0.369

#subset y 
ind<-which(rownames(hmscYlo3)%in%rownames(hmscXlo2))
hmscYlo4<-hmscYlo3[ind,]

#subset dist
ind<-which(comm.bio$X.SampleID%in%rownames(hmscXlo2))
hmiscDISTlo<-data.frame(X=comm.bio$X,Y=comm.bio$Y,elevation=comm.bio$elevation)[ind,]
rownames(hmiscDISTlo)<-comm.bio$X.SampleID[ind]

#sort them the same
hmscYlo5<-hmscYlo4[order(rownames(hmscYlo4)),]
hmscXlo3<-hmscXlo2[order(rownames(hmscXlo2)),]
hmiscDISTlo2<-hmiscDISTlo[order(rownames(hmiscDISTlo)),]


#me
ind<-which(comm.bio$lomehi=="me")
hmscXme<-data.frame(snowdepth=comm.bio$snowdepth,pH=comm.bio$pH,moisture=comm.bio$moisture,cvsnow=comm.bio$cvsnow,pca1=comm.bio$pca1)[ind,]
rownames(hmscXme)<-comm.bio$X.SampleID[ind]

which(hmscXme$pca1<(-0.330213)|(hmscXme$pca1>(-0.012)))
set.seed(6);ind<-c(sample(c(1:6,8:16,18:25),size=10),7,17);ind #make sure 7 and 17 are not included
hmscXme2<-data.frame(snowdepth=hmscXme$snowdepth,pH=hmscXme$pH,moisture=hmscXme$moisture,cvsnow=hmscXme$cvsnow)[ind,]
rownames(hmscXme2)<-rownames(hmscXme)[ind]
#range is 0.3192

#subset y 
ind<-which(rownames(hmscYme3)%in%rownames(hmscXme2))
hmscYme4<-hmscYme3[ind,]

#subset dist
ind<-which(comm.bio$X.SampleID%in%rownames(hmscXme2))
hmiscDISTme<-data.frame(X=comm.bio$X,Y=comm.bio$Y,elevation=comm.bio$elevation)[ind,]
rownames(hmiscDISTme)<-comm.bio$X.SampleID[ind]

#sort them the same
hmscYme5<-hmscYme4[order(rownames(hmscYme4)),]
hmscXme3<-hmscXme2[order(rownames(hmscXme2)),]
hmiscDISTme2<-hmiscDISTme[order(rownames(hmiscDISTme)),]


#hi
ind<-which(comm.bio$lomehi=="hi")
hmscXhi<-data.frame(snowdepth=comm.bio$snowdepth,pH=comm.bio$pH,moisture=comm.bio$moisture,cvsnow=comm.bio$cvsnow,pca1=comm.bio$pca1)[ind,]
rownames(hmscXhi)<-comm.bio$X.SampleID[ind]
ind<-which(hmscXhi$pca1>0.1074467&hmscXhi$pca1<0.4914467)
hmscXhi2<-hmscXhi[ind,1:4]
#range is 0.383

#subset y
ind<-which(rownames(hmscYhi3)%in%rownames(hmscXhi2))
hmscYhi4<-hmscYhi3[ind,]

#subset dist
ind<-which(comm.bio$X.SampleID%in%rownames(hmscXhi2))
hmiscDISThi<-data.frame(X=comm.bio$X,Y=comm.bio$Y,elevation=comm.bio$elevation)[ind,]
rownames(hmiscDISThi)<-comm.bio$X.SampleID[ind]

#sort them the same, should be the same, but just to be crazy
hmscYhi5<-hmscYhi4[order(rownames(hmscYhi4)),]
hmscXhi3<-hmscXhi2[order(rownames(hmscXhi2)),]
hmiscDISThi2<-hmiscDISThi[order(rownames(hmiscDISThi)),]


##### Select common species #####
#select species with greater than X (X+1 or more) occurrences

## low ##
ind<-which(colSums(hmscYlo5>0)>8)
#hmscYlo2sub<-hmscYlo2[,ind] #subset the raw data as well just for looking at histograms
length(ind)
hmscYlo6<-hmscYlo5[,ind]
dim(hmscYlo6)

## medium ##
ind<-which(colSums(hmscYme5>0)>8)
length(ind)
hmscYme6<-hmscYme5[,ind]
hmscXme
dim(hmscYme6)

## high ##
ind<-which(colSums(hmscYhi5>0)>8)
length(ind)
hmscYhi6<-hmscYhi5[,ind]
dim(hmscYhi6)

hmscXlo2<-scale(hmscXlo)
hmscXme2<-scale(hmscXme)
hmscXhi2<-scale(hmscXhi)

#Make them matrices
hmscYlo7<-as.matrix(hmscYlo6)
hmscYme7<-as.matrix(hmscYme6)
hmscYhi7<-as.matrix(hmscYhi6)

hmscXlo4<-as.matrix(hmscXlo3)
hmscXme4<-as.matrix(hmscXme3)
hmscXhi4<-as.matrix(hmscXhi3)

hmiscDISTlom<-as.matrix(dist(hmiscDISTlo2))
hmiscDISTmem<-as.matrix(dist(hmiscDISTme))
hmiscDISThim<-as.matrix(dist(hmiscDISThi2))


#### Modeling ####
dim(hmscYlo7)
mod.lo<- boral(y = hmscYlo7, X = hmscXlo4, lv.control = list(num.lv = 2,type="spherical",distmat=hmiscDISTlom), family = c(rep("normal",106),rep("negative.binomial",0)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.lo <- get.residual.cor(mod.lo) 

dim(hmscYme7)
mod.me<- boral(y = hmscYme7, X = hmscXme4, lv.control = list(num.lv = 2,type="spherical",distmat=hmiscDISTmem), family = c(rep("normal",89),rep("negative.binomial",1)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.me <- get.residual.cor(mod.me) 

dim(hmscYhi7)
mod.hi<- boral(y = hmscYhi7, X = hmscXhi4, lv.control = list(num.lv = 2,type="spherical",distmat=hmiscDISThim), family = c(rep("normal",96),rep("negative.binomial",3)), save.model = TRUE, calc.ics = T, mcmc.control = list(n.burnin = 10000, n.iteration = 40000, n.thin = 30, seed = 123))#
rescor.hi <- get.residual.cor(mod.hi) 



# >8, 3lv, complexity good, total interactions not
# >8, 2lv, complexity good, total interactions good, let's use this






##### Network diagrams for lo #####
#creating sparse matrix
colMatlo<-rescor.lo$sig.correlaton#rescor.lo11auto
colMatlo[which(rescor.lo$sig.correlaton>0)]<-1
colMatlo[which(rescor.lo$sig.correlaton<0)]<- -1

graphlo1<-graph_from_adjacency_matrix(colMatlo, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistlo<-data.frame(as_edgelist(graphlo1),weight=E(graphlo1)$weight) #just the edges

length(which(myedgelistlo$weight==1))
length(which(myedgelistlo$weight==-1))
length(which(myedgelistlo$weight==1))/(length(which(myedgelistlo$weight==1))+length(which(myedgelistlo$weight==-1)))

graphlo2<-graph.edgelist(as.matrix(myedgelistlo[,1:2]),directed=FALSE)
graphlo2

verticesgraphlo<-data.frame(otu=rownames(as.matrix(V(graphlo2))))
colorgraphlo<-merge(verticesgraphlo,labelsall,"otu",all.y=F,all.x=F,sort=F)
#shapesgraplo<-ifelse(colorgraphlo$group%in%c("Eukaryota"),"csquare",'circle')

##use colorgraphlo$group2 for ordering, if ther are ties it leaves them in thir original order, thus is still preserves some of the ordering that makes the lines look nice
orderlo<-order(colorgraphlo$group2)
#orderlo<-order(verticesgraphlo$otu)
graphlo2$layout <- layout_in_circle(graphlo2,order=orderlo)
#graphlo2$layout <- layout_in_circle

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/networklocircleposblue.pdf") 
#plot(graphlo2,vertex.size=4,edge.curved=F,edge.color=ifelse(myedgelistlo$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphlo$color,edge.width=.7,vertex.label=NA)#,vertex.shape=shapesgraplo  positive is red
plot(graphlo2,vertex.size=4,edge.curved=F,edge.color=ifelse(myedgelistlo$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphlo$color,edge.width=.7,vertex.label=NA)#,vertex.shape=shapesgraplo  positive is blue
#dev.off()





##### Network diagrams for me #####
#creating sparse matrix
colMatme<-rescor.me$sig.correlaton
colMatme[which(rescor.me$sig.correlaton>0)]<-1
colMatme[which(rescor.me$sig.correlaton<0)]<- -1

graphme1<-graph_from_adjacency_matrix(colMatme, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelistme<-data.frame(as_edgelist(graphme1),weight=E(graphme1)$weight) #just the edges

length(which(myedgelistme$weight==1))
length(which(myedgelistme$weight==-1))
length(which(myedgelistme$weight==1))/(length(which(myedgelistme$weight==1))+length(which(myedgelistme$weight==-1)))

graphme2<-graph.edgelist(as.matrix(myedgelistme[,1:2]),directed=FALSE)
graphme2

verticesgraphme<-data.frame(otu=rownames(as.matrix(V(graphme2))))
colorgraphme<-merge(verticesgraphme,labelsall,"otu",all.y=F,all.x=F,sort=F)

##use colorgraphlo$group2 for ordering, if there are ties it leaves them in their original order, thus is still preserves some of the ordering that makes the lines look nice
orderme<-order(colorgraphme$group2)
#orderme<-order(verticesgraphme$otu)
graphme2$layout <- layout_in_circle(graphme2,order=orderme)
#graphme2$layout <- layout_in_circle

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/networkmecircleposblue.pdf") 
#plot(graphme2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistme$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphme$color,edge.width=.7)#,layout=l3
plot(graphme2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelistme$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphme$color,edge.width=.7)#,layout=l3
#dev.off()



##### Network diagrams for hi #####
#creating sparse matrix
colMathi<-rescor.hi$sig.correlaton#rescor.hi11flv3
colMathi[which(rescor.hi$sig.correlaton>0)]<-1
colMathi[which(rescor.hi$sig.correlaton<0)]<- -1

graphhi1<-graph_from_adjacency_matrix(colMathi, mode = c( "undirected"), weighted = T, diag = F,add.colnames = NULL, add.rownames = NULL)
myedgelisthi<-data.frame(as_edgelist(graphhi1),weight=E(graphhi1)$weight) #just the edges

length(which(myedgelisthi$weight==1))
length(which(myedgelisthi$weight==-1))
length(which(myedgelisthi$weight==1))/(length(which(myedgelisthi$weight==1))+length(which(myedgelisthi$weight==-1)))

graphhi2<-graph.edgelist(as.matrix(myedgelisthi[,1:2]),directed=FALSE)
graphhi2

verticesgraphhi<-data.frame(otu=rownames(as.matrix(V(graphhi2))))
colorgraphhi<-merge(verticesgraphhi,labelsall,"otu",all.y=F,all.x=F,sort=F)

#order starts at 3:00 and goes counterclockwise
##use colorgraphlo$group2 for ordering, if ther are ties it leaves them in their original order, thus is still preserves some of the ordering that makes the lines look nice
orderhi<-order(colorgraphhi$group2)
#orderhi<-order(verticesgraphhi$otu)
graphhi2$layout <- layout_in_circle(graphhi2,order=orderhi)
#graphhi2$layout <- layout_in_circle

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/networkhicircleposblue.pdf") 
#plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#ce4d42","#687dcb"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
plot(graphhi2,vertex.size=4,edge.curved=F,vertex.label=NA,edge.color=ifelse(myedgelisthi$weight==1,"#687dcb","#ce4d42"),vertex.color=colorgraphhi$color,edge.width=.7)#,layout=l3  
#dev.off()



length(E(graphlo2))/length(V(graphlo2))
length(E(graphme2))/length(V(graphme2))
length(E(graphhi2))/length(V(graphhi2))

length(E(graphlo2))
length(E(graphme2))
length(E(graphhi2))
