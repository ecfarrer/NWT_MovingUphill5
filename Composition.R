###### Change in relative abundance ######
datBacS5k2
datITSS5k2
datEukS5k2
datEukN5k2

#change two colnames, they can't have a dash in it
colnames(datBacS5k2)[85]<-"WPS.2"
colnames(datBacS5k2)[38]<-"BHI80.139"
names(which(colSums(datBacS5k2[,34:88])>2))
relBac<-datBacS5k2 %>% 
  dplyr::select(Sample_name,Acidobacteria,Bacteroidetes,Cyanobacteria,Gemmatimonadetes,Heterotrophic_Actinobacteria,Heterotrophic_Chloroflexi,Heterotrophic_Planctomycetes,Heterotrophic_Proteobacteria,Heterotrophic_Verrucomicrobia) %>% #,WPS.2 is unknown function so I'm taking it out
  gather(Taxa,abun,Acidobacteria:Heterotrophic_Verrucomicrobia) %>%
  mutate(Taxa = recode_factor(Taxa, Acidobacteria="Acidobacteria (H)",Bacteroidetes = "Bacteroidetes (H)",Cyanobacteria="Cyanobacteria (P)",Gemmatimonadetes="Gemmatimonadetes (H)",Heterotrophic_Actinobacteria="Actinobacteria (H)",Heterotrophic_Chloroflexi="Chloroflexi (H)",Heterotrophic_Planctomycetes="Planctomycetes (H)",Heterotrophic_Proteobacteria="Proteobacteria (H)",Heterotrophic_Verrucomicrobia="Verrucomicrobia (H)"))%>%
  mutate(type="A. Bacteria")

names(which(colSums(datITSS5k2[,34:47])>.75))
relITS<-datITSS5k2 %>% 
  dplyr::select(Sample_name,Ascomycota,Basidiomycota,Glomeromycota,Mortierellomycota) %>%
  gather(Taxa,abun,Ascomycota:Mortierellomycota) %>%
  mutate(type="B. Fungi")

#if the label has the word "unknown" in it, then I don't want to plot it. it means that is it unknown at a level higher than phylum (even if I could tell hetero/photo)
sort(colSums(datEukS5k2[,34:62]))
names(which(colSums(datEukS5k2[,34:62])>1))#,Alveolata,Archaeplastida,Photosynthetic_Stramenopiles,Rhizaria
relEukS<-datEukS5k2 %>% 
  dplyr::select(Sample_name,Cercozoa,Charophyta,Chlorophyta,Ciliophora,Heterotrophic_Euglenozoa,Heterotrophic_Stramenopiles,Photosynthetic_Stramenopiles) %>%
  gather(Taxa,abun,Cercozoa,Charophyta,Chlorophyta,Ciliophora,Heterotrophic_Euglenozoa,Heterotrophic_Stramenopiles,Photosynthetic_Stramenopiles) %>%
  mutate(Taxa = recode_factor(Taxa, Cercozoa="Cercozoa (H)",Charophyta = "Charophyta (P)",Chlorophyta="Chlorophyta (P)",Ciliophora="Ciliophora (H)",Heterotrophic_Euglenozoa="Euglenozoa (H)",Heterotrophic_Stramenopiles="Stramenopiles (H)",Photosynthetic_Stramenopiles="Stramenopiles (P)"))%>%
  mutate(type="C. Small Eukaryotes")

names(which(colSums(datEukN5k2[,34:41])>1))
relEukN<-datEukN5k2 %>% 
  dplyr::select(Sample_name,Arthropoda,Nematoda,Rotifera,Tardigrada) %>%
  gather(Taxa,abun,Arthropoda,Nematoda,Rotifera,Tardigrada) %>%
  mutate(type="D. Soil microfauna")

relALL1<-rbind(relBac,relITS,relEukS,relEukN)#
head(relALL1)

#merge with biogeo6 to get pca1
relALL<-merge(relALL1,biogeo6,"Sample_name")
head(relALL)

#plotdata<-relALL %>%
#  mutate(typeTaxa=paste(type,Taxa)) %>%
#  group_by(Taxa,lomehi,type,typeTaxa) %>%
#  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
#  #%>%filter(mean_abun>.04)

#this was weird, maybe something changed in ggplot or dplyr because the colors were messing up and it was listing the legend in alfabetical order by taxa rather than the order in the "plotdata" dataframe. the workaroudn was to set the levels of plotdata$Taxa so they were correct
plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa)) %>%
  merge(anovaoutputauto2)%>%
  group_by(typeTaxa,Taxa,lomehi,type,sig) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun))
plotdata$Taxa<-factor(plotdata$Taxa,levels=unique(plotdata$Taxa))

mylevels<-levels(plotdata$Taxa)

as.data.frame(plotdata)
plotdata$lomehi<-factor(plotdata$lomehi,levels=c("lo","me","hi"))

#9 bacteria, 4 fungi, 7 small euks, 4 large euks
mycols<-c("#D9A125",#yellow
          #"#4BC366",#light green
          "#6F94DE",#light blue
          "#B4405E",#red
          "#D185E0",#light purple
          "#659125",#green
          "#ff99a4",#light pink
          "#cf6f23",#orange
          "#5C426C",#dark purple
          "#6768A3",#medium blue last bact
          
          "#cf6f23",#orange" 
          "#D9A125",#yellow
          "#B4405E",#red
          "#6768A3",#medium blue
          
          "#cf6f23",#orange
          "#659125",#green
          "#4BC366",#light green          
          "#D185E0",#light purple
          "#D9A125",#yellow
          "#5C426C",#dark purple 
          "#008B8B",#greenblue "#ff99a4",#light pink       
          
          "#cf6f23",#orange
          "#D9A125", #yellow          
          "#5C426C", #dark purple
          "#6F94DE")#light blue

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/relabuntaxavsplantdensitygroupsR2.pdf",width=6.5,height=4.3)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=lomehi,y=mean_abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5, aes(linetype=sig))+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.15,size=.5)+
  scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
#dev.off()



#scatter plots - super messy
head(relALL)

#need to run the linear models below first to get lmoutput2
plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa))%>%
  merge(lmoutput2)
plotdata$Taxa<-factor(plotdata$Taxa,levels=mylevels)


#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforFrontiersSubmission/relabuntaxavsplantdensitygroupsBFSLENscatter.pdf",width=6.5,height=4.3)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=log(pca1+1),y=abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  scale_color_manual(values=mycols) +
  geom_point(size=.1,shape=16)+
  #geom_smooth(method=lm,se=F,size=.8) +
  geom_smooth(method = "lm", formula = y ~ poly(x,2),se=F,size=.5, aes(linetype=sig)) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()


#Doing anova on all of the above taxa groups
ind<-length(unique(relALL$Taxa))
anovaoutput<-data.frame(Taxa=rep(NA,ind),F=rep(NA,ind),P=rep(NA,ind))
for(i in 1:ind){
  current.taxa<-levels(relALL$Taxa)[i]
  temp<-relALL %>%
    filter(Taxa==current.taxa)
  mod<-anova(lm(abun~lomehi,data=temp))
  anovaoutput[i,1]<-as.character(current.taxa)
  anovaoutput[i,2]<-mod$`F value`[1]
  anovaoutput[i,3]<-mod$`Pr(>F)`[1]
}
anovaoutput$qval<-p.adjust(anovaoutput$P,method="fdr")
anovaoutput$qval<-format(anovaoutput$qval,scientific=F)

#anovaoutput$Taxa<-factor(anovaoutput$Taxa,levels=unique(plotdata$Taxa))
anovaoutput[order(anovaoutput$Taxa),]


#Doing anova including spatial autocorrelation on all of the above taxa groups
ind<-length(unique(relALL$Taxa))
anovaoutputauto<-data.frame(Taxa=rep(NA,ind),F=rep(NA,ind),P=rep(NA,ind))
for(i in 1:ind){
  current.taxa<-levels(relALL$Taxa)[i]
  temp<-relALL %>%
    filter(Taxa==current.taxa)
  mod2<-gls(abun~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=temp)
  amod2<-anova(mod2)
  anovaoutputauto[i,1]<-as.character(current.taxa)
  anovaoutputauto[i,2]<-amod2$`F-value`[2]
  anovaoutputauto[i,3]<-amod2$`p-value`[2]
}
anovaoutputauto$qval<-p.adjust(anovaoutputauto$P,method="fdr")
anovaoutputauto$qval<-format(anovaoutputauto$qval,scientific=F)
anovaoutputauto$sig<-as.factor(ifelse(anovaoutputauto$qval<0.1,"sig","signot"))
anovaoutputauto2<-anovaoutputauto%>%select(Taxa,sig)





#Doing (non)linear regression on all of the above taxa groups
ind<-length(unique(relALL$Taxa))
lmoutput<-data.frame(Taxa=rep(NA,ind),P=rep(NA,ind))
relALL$lpca1<-log(relALL$pca1+1)
relALL$lpca12<-(log(relALL$pca1+1))^2
for(i in 1:ind){
  current.taxa<-levels(relALL$Taxa)[i]
  temp<-relALL %>%
    filter(Taxa==current.taxa)
  mod2<-gls(abun~lpca1+lpca12,correlation=corSpher(form = ~ X + Y + elevation),data=temp,method="ML")
  mod1<-gls(abun~1,correlation=corSpher(form = ~ X + Y + elevation),data=temp,method="ML")
  compmod<-anova(mod1,mod2)
  lmoutput[i,1]<-as.character(current.taxa)
  lmoutput[i,2]<-compmod$`p-value`[2]
}
lmoutput$qval<-p.adjust(lmoutput$P,method="fdr")
lmoutput$qval<-format(lmoutput$qval,scientific=F)
lmoutput$sig<-as.factor(ifelse(lmoutput$qval<0.05,"sig","signot"))
#lmoutput$linetype<-ifelse(lmoutput$P<0.05,"sig","nonsig")
lmoutput2<-lmoutput%>%select(Taxa,sig)

plot(temp$lpca1,temp$abun)
curve(.16+0.024*x+0.00517*x,add=T)

#anovaoutput$Taxa<-factor(anovaoutput$Taxa,levels=unique(plotdata$Taxa))
anovaoutput[order(anovaoutput$Taxa),]






##### looking at patterns of rarer groups #####
relBacsub<-datBacS5k2 %>% 
  dplyr::select(Sample_name,AD3:WS5)%>%
  gather(Taxa,abun,AD3:WS5) 
relBacsub2<-merge(relBacsub,biogeo6,"Sample_name")
head(relBacsub2)
plotdatasub<-relBacsub2 %>%
  group_by(Taxa,lomehi) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
plotdatasub$Taxa<-factor(plotdatasub$Taxa,levels=unique(plotdatasub$Taxa))
plotdatasub$lomehi<-factor(plotdatasub$lomehi,levels=c("lo","me","hi"))

data.frame(plotdatasub)






##### Ordination #####

#input file, from beginning of boral script
#comm.bio is nice b/c it is already merged with biogeo6 so the microbial data have been reduced to 75 total samples (vs. 90) and b/c all microbial datasets are then all in the same order
comm.ord<-comm.bio

#calculate CLR on full datasets for 16S, ITS, EukS, EukN

colnames(comm.ord)
temp<-lapply(colnames(comm.ord),function(x){substr(x,1,1)})
head(which(temp=="I"))
tail(which(temp=="I"))

#separate N, S, 16S, ITS
env<-comm.ord[,1:54]
commN<-comm.ord[,55:504]
commS<-comm.ord[,505:3665]
commB<-comm.ord[,3666:20277]
commI<-comm.ord[,20278:24511]

rownames(commN)<-env$X.SampleID
rownames(commS)<-env$X.SampleID
rownames(commB)<-env$X.SampleID
rownames(commI)<-env$X.SampleID

#take out doubletons and singletons
ind<-which(colSums(commN>0)>2);length(ind)
commN2<-commN[ind]
ind<-which(colSums(commS>0)>2);length(ind)
commS2<-commS[ind]
ind<-which(colSums(commB>0)>2);length(ind)
commB2<-commB[ind]
ind<-which(colSums(commI>0)>2);length(ind)
commI2<-commI[ind]

#estimate zeros and calculate clr
commN3 <- cmultRepl(commN2,label=0, method="CZM")
commN4 <- t(apply(commN3, 1, function(x){log(x) - mean(log(x))}))
commS3 <- cmultRepl(commS2,label=0, method="CZM")
commS4 <- t(apply(commS3, 1, function(x){log(x) - mean(log(x))}))
commB3 <- cmultRepl(commB2,label=0, method="CZM")
commB4 <- t(apply(commB3, 1, function(x){log(x) - mean(log(x))}))
commI3 <- cmultRepl(commI2,label=0, method="CZM")
commI4 <- t(apply(commI3, 1, function(x){log(x) - mean(log(x))}))

#just the stats:
adonis2(commB4~lomehi,data=env,permutations = how(nperm=10000),method="euclidean")
adonis2(commI4~lomehi,data=env,permutations = how(nperm=10000),method="euclidean")
adonis2(commS4~lomehi,data=env,permutations = how(nperm=10000),method="euclidean")
adonis2(commN4~lomehi,data=env,permutations = how(nperm=10000),method="euclidean")

#plots and extra analyses and visuals are below

#pca
ordN <- prcomp(commN4)
totvarN<-sum(ordN$sdev^2)
PC1 <- paste("PC1: ", round(sum(ordN$sdev[1]^2)/totvarN, 3))
PC2 <- paste("PC2: ", round(sum(ordN$sdev[2]^2)/totvarN, 3))
biplot(ordN, var.axes=T, scale=0, xlab=PC1, ylab=PC2)
plot(scores(ordN),col=col,bg=col,pch=21,cex=2)

#the gloor paper suggests using pca then anosim
anosim(commN4,grouping=env$lomehi,distance="euclidean") #pretty much same as adonis, but adonis is more robust. adonis I think is the same as anova but maybe with slightly different math
adonis2(commN4~lomehi,data=env,permutations = how(nperm=10000),method="euclidean")

#rda
myrda <- rda(commN4 ~ lomehi, data = env)
summary(myrda)
anova(myrda, by="margin", permutations = 1000)
plot(myrda,scaling=1)


#### dbrda - dbrda is same as rda with euclidean distance ####

#nematodes, explains 8.2% p=.001
mydbrda<-dbrda(commN4~lomehi,distance="euclidean",data=env)
anova(mydbrda, by = "margin",permutations = how(nperm=10000))

col<-ifelse(env$lomehi=="lo","#9350a1",ifelse(env$lomehi=="me","#697cd4",ifelse(env$lomehi=="hi","#62ad64","#b8475f")))
#low purple, meblue , hi green, else red

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/ordinationN.pdf",width=3.8,height=4.2)
plot(scores(mydbrda)$sites,col=col,bg=col,pch=21,cex=1,cex.lab=.8,cex.axis=.8,xlab="Axis 1 (5.4%)",ylab="Axis 2 (2.7%)")#,xlim=c(-4,3.5),ylim=c(-5,4)
#text(scores(mydbrda)$centroids,labels=c("Late","Early","Mid"),col=c("#62ad64","#9350a1","#697cd4"),cex=2)
#text(scores(mydbrda)$sites,labels=env$Sample_name)
#the ellipses look crappy and big
#ordiellipse(mydbrda,groups=col,col=c("#62ad64","#697cd4","#9350a1"),conf=.95,kind="sd",lwd=2)#
legend("bottomleft",c("Early","Mid","Late"),cex=.8,pch=21,col=c("#9350a1","#697cd4","#62ad64"),pt.bg=c("#9350a1","#697cd4","#62ad64"),bty="n")
dev.off()


#Small euks, explains 8.2% p=0.001
mydbrda<-dbrda(commS4~lomehi,distance="euclidean",data=env)
anova(mydbrda, by = "margin",permutations = how(nperm=10000))

col<-ifelse(env$lomehi=="lo","#9350a1",ifelse(env$lomehi=="me","#697cd4",ifelse(env$lomehi=="hi","#62ad64","#b8475f")))
#low purple, meblue , hi green, else red

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/ordinationS.pdf",width=3.8,height=4.2)
plot(scores(mydbrda)$sites,col=col,bg=col,pch=21,cex=1,cex.lab=.8,cex.axis=.8,xlab="Axis 1 (6.4%)",ylab="Axis 2 (1.8%)")
#text(scores(mydbrda)$centroids,labels=c("High","Low","Mid"),col=c("#62ad64","#9350a1","#697cd4"),cex=2)
#text(scores(mydbrda)$sites,labels=env$Sample_name)
legend("bottomleft",c("Early","Mid","Late"),cex=.8,pch=21,col=c("#9350a1","#697cd4","#62ad64"),pt.bg=c("#9350a1","#697cd4","#62ad64"),bty="n")
dev.off()


#Bact, explains 9.3%, p=0.001
mydbrda<-dbrda(commB4~lomehi,distance="euclidean",data=env)
anova(mydbrda, by = "margin",permutations = how(nperm=10000))

col<-ifelse(env$lomehi=="lo","#9350a1",ifelse(env$lomehi=="me","#697cd4",ifelse(env$lomehi=="hi","#62ad64","#b8475f")))
#low purple, meblue , hi green, else red

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/ordinationB.pdf",width=3.8,height=4.2)
plot(scores(mydbrda)$sites,col=col,bg=col,pch=21,cex=1,cex.lab=.8,cex.axis=.8,xlab="Axis 1 (7.4%)",ylab="Axis 2 (2.0%)")
#text(scores(mydbrda)$centroids,labels=c("High","Low","Mid"),col=c("#62ad64","#9350a1","#697cd4"),cex=2)
#text(scores(mydbrda)$sites,labels=env$Sample_name)
legend("bottomleft",c("Early","Mid","Late"),cex=.8,pch=21,col=c("#9350a1","#697cd4","#62ad64"),pt.bg=c("#9350a1","#697cd4","#62ad64"),bty="n")
dev.off()

#Fungi, explains 10.2% p=0.001
mydbrda<-dbrda(commI4~lomehi,distance="euclidean",data=env)
anova(mydbrda, by = "margin",permutations = how(nperm=10000))

col<-ifelse(env$lomehi=="lo","#9350a1",ifelse(env$lomehi=="me","#697cd4",ifelse(env$lomehi=="hi","#62ad64","#b8475f")))
#low purple, meblue , hi green, else red

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/ordinationI.pdf",width=3.8,height=4.2)
plot(scores(mydbrda)$sites,col=col,bg=col,pch=21,cex=1,cex.lab=.8,cex.axis=.8,xlab="Axis 1 (8.3%)",ylab="Axis 2 (1.9%)")
#text(scores(mydbrda)$centroids,labels=c("High","Low","Mid"),col=c("#62ad64","#9350a1","#697cd4"),cex=2)
#text(scores(mydbrda)$sites,labels=env$Sample_name)
legend("bottomright",c("Early","Mid","Late"),cex=.8,pch=21,col=c("#9350a1","#697cd4","#62ad64"),pt.bg=c("#9350a1","#697cd4","#62ad64"),bty="n")
dev.off()











