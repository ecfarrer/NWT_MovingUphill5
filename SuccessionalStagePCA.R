#Successional stage PCA


##### biogeochemistry #####
biogeo<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Biogeochemistry/CN_enzymes_pH_moisture_whc_Niwot2015t.csv")
head(biogeo)
dim(biogeo)
#note on dataset - I updated this because for the data from microbial biomass (gravimetric moisture, IN, DOC, etc), sample 33 and 34 were mixed up. origianlly sample 33 was missing from the first dataset dorota gave me (CN_enzymes_pH_moisture_whc_Niwot2015.xlsx), and the nubmers fom 34 were really the values for 33. Then I got the All Data file from Dorota and this cleared it up b/c it had both samples 33 adn 34 in it, I also double checked with the Biogeochemistry_Niwot_2015.xlsx file where the calculations were done.

#merge with one of the mapping files
biogeo2<-merge(biogeo,datEukS5otu[,1:31])
#cbind(biogeo2$moisture,biogeo2$WHC)
head(biogeo2)

#I could replace all negative numbers with 0, not sure if I should do this, since they are relative, but it doesn't really make sense to have negative microbial biomass or inorganic N
#biogeo2[biogeo2<0]<-0


#king snowdepth (to get cv and 2015 snowdepth)
snowdepth<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/King_SnowDepth.csv")
head(snowdepth)
colnames(snowdepth)[1]<-"Sample_name"
colnames(snowdepth)[24]<-"snow2015"
snowdepth$cvsnow<-apply(snowdepth[,6:24],1,function(x){sd(x,na.rm=T)/mean(x,na.rm=T)})
snowdepthr<-dplyr::select(snowdepth,Sample_name,snow2015,cvsnow)

biogeo3<-merge(biogeo2,snowdepthr,"Sample_name") #with merge, they got put in the right order
rownames(biogeo3)<-biogeo3$X.SampleID
head(biogeo3)



##### plant cover #####
#get plant cover data to use as "light" and merge with biogeo
plantcov<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Plants/Niwot_MovingUpHill_plots2015.csv")
names(plantcov)[1]<-"Sample_name"
plantcov$plantcover<-(plantcov$MOSS+plantcov$VEG)/plantcov$TOTAL
plantcov$plantcover[which(plantcov$Sample_name==64)]<-0 #this I'm 100% sure has zero plants, however bare and rock were recorded as 0s. I looked up the sample data sheet and it looks like Sam and I did not sample it because a soil smaple was not taken in 2007, however since it had 0 plants we added that plot and forgot to do cover on it.
head(plantcov)
plantcov$vascularplantcover<-(plantcov$VEG)/plantcov$TOTAL
plantcov$vascularplantcover[which(plantcov$Sample_name==64)]<-0

#take out some columns that I don't want to confuse with mapping file columns
plantcov2<-plantcov%>%select(Sample_name,X,Y,plantcover,vascularplantcover)

head(biogeo3)
biogeo4<-merge(biogeo3,plantcov2,"Sample_name") #with merge, they got put in the right order
rownames(biogeo4)<-biogeo4$X.SampleID
head(biogeo4)
#write.csv(biogeo3,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/biogeo3.csv",row.names=F)



##### PCA #####
#Use these plant variables: they include mosses/liverworts but not lichen
#Plant_Dens
#Plant_Div

#datITSS3otu3$X

biogeo4a<-biogeo4%>%
  select(Plant_Dens,Plant_Div,TC,TN,NH4,NO3,pH,WHC,moisture,snowdepth,elevation,plantcover,MicC,MicN)#,IN,DOC,DON,
ind<-which(is.na(rowSums(biogeo4a))==F)
biogeo5<-biogeo4a[ind,]
mynames<-rownames(biogeo5)

# #for correlation matrix
# biogeo4a<-biogeo4%>%
#  select(TC,TN,NH4,NO3,MicC,MicN,pH,WHC,moisture,snowdepth,elevation,Plant_Dens,Plant_Div,plantcover,snow2015,cvsnow)#,IN,DOC,DON,
# ind<-which(is.na(rowSums(biogeo4a))==F)
# ind<-which(rownames(biogeo4a)%in%mynames)
# biogeo5<-biogeo4a[ind,]
# dim(biogeo5)
# write.csv(round(rcorr(as.matrix(biogeo5))$r,2),"corrR.csv")
# write.csv(round(rcorr(as.matrix(biogeo5))$P,4),"corrP.csv")

#trying to log transform all skewed variables:
#pH, snowdepth,elevation is more or less normal
# biogeo5$Plant_Dens<-log(biogeo5$Plant_Dens+1)
# biogeo5$Plant_Div<-log(biogeo5$Plant_Div+1)
# biogeo5$TC<-log(biogeo5$TC)
# biogeo5$TN<-log(biogeo5$TN)
# biogeo5$NH4<-log(biogeo5$NH4+.2)
# biogeo5$NO3<-log(biogeo5$NO3+.2)
# biogeo5$WHC<-log(biogeo5$WHC)
# biogeo5$moisture<-log(biogeo5$moisture)
# biogeo5$plantcover<-log(biogeo5$plantcover+.02)
# biogeo5$MicC<-log(biogeo5$MicC+150)
# biogeo5$MicN<-log(biogeo5$MicN+10)
#surprisingly this doesn't help too much, the range of lo is 0.83, me 0.56, hi 1.8 (two fold difference)
#the regular way the range of lo is 0.38, me 0.32, hi 2.9 (7 fold difference)

mypca<-rda(biogeo5,scale=T,na.action=na.omit)
#plot(mypca)
#summary(mypca)

scal=2
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/successionpca.pdf")
plot(0,0,type="n",xlab=paste("Axis 1 (",sprintf("%.1f",mypca$CA$eig["PC1"]/mypca$tot.chi*100,3),"%)",sep=""),ylab=paste("Axis 2 (",sprintf("%.1f",mypca$CA$eig["PC2"]/mypca$tot.chi*100,3),"%)",sep=""),cex.lab=1.4,xlim=c(-1,3),ylim=c(-1.5,1.5))
abline(h=0)
abline(v=0)
sorts<-order(scores(mypca,scaling=scal)$sites[,1],decreasing=T)
points(-scores(mypca,scaling=scal)$sites[sorts,1],scores(mypca,scaling=scal)$sites[sorts,2],cex=.9,pch=16,col=rep(c('gray10','gray50','gray80'),each=25))
arrows(0,0,-scores(mypca,scaling=scal)$species[,1],scores(mypca,scaling=scal)$species[,2],length=.05,col=1,lwd=2)
text(-scores(mypca,scaling=scal)$species[,1],scores(mypca,scaling=scal)$species[,2],labels=rownames(scores(mypca,scaling=scal)$species),cex=.9,col=1)
legend(2,1.5,c("Early succession","Mid succession","Late succession"),pch=16,col=c('gray10','gray50','gray80'),bty="n",cex=.9)
#dev.off()


succession<-(-scores(mypca)$sites[,1])

succession2<-data.frame(pca1=succession,X.SampleID=names(succession))
succession3<-succession2[order(succession2$pca1),]
#succession3old<-succession2[order(succession2$pca1),]
succession3$lomehi<-rep(c('lo','me','hi'),each=25)
#succession3old$lomehi<-rep(c('lo','me','hi'),each=25)
hist(succession3$pca1)

#cbind(succession3,succession3old) #when I log transform, the plot classification into lo me hi is different

plot(biogeo2$snowdepth,biogeo2$Plant_Dens)
plot(biogeo2$snowdepth,biogeo2$pH)
plot(biogeo2$Plant_Dens,biogeo2$pH)

summary(lm(biogeo3$Plant_Dens~biogeo3$WHC))
summary(lm(biogeo3$Plant_Dens~biogeo3$moisture))


###### Merge back with biogeo3 #####
biogeo6<-merge(biogeo4[,-which(names(biogeo4)=="lomehi")],succession3)
head(biogeo6)
cbind(biogeo6$plantcover,biogeo6$Plant_Dens,biogeo6$lomehi)
#delete two columns so it is the same dimensions as before but with snow2015 and cvsnow
biogeo6$Description<-NULL
biogeo6$BarcodeSequence<-NULL
#it is kind of strange that some of the plots with 0 plants were classified as medium, but there is not much I can do if that's how the pca shakes out. i'll just have to see what the networks look like


print(rcorr(as.matrix(biogeo5))$r,digits=10)
print(rcorr(biogeo6$snowdepth,biogeo6$cvsnow)$r,digits=10)

#write file for Cliff to make a map of the plots
#write.csv(biogeo6,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/networkmsplots.csv")





