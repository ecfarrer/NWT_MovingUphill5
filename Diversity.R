
##### Phylogenetic Diversity #####
#Files:
richEukS2
richEukN2
richBac2
richITS2

biogeo6

richBac2$X.SampleID<-rownames(richBac2)
richBac3<-merge(richBac2,biogeo6,"X.SampleID")
richBac3$type<-"1Bacteria"

richITS2$X.SampleID<-rownames(richITS2)
richITS3<-merge(richITS2,biogeo6,"X.SampleID")
richITS3$type<-"2Fungi"

#use chao1 for fungi
richITS3$PD<-richITS3$Chao1

richEukS2$X.SampleID<-rownames(richEukS2)
richEukS3<-merge(richEukS2,biogeo6,"X.SampleID")
richEukS3$type<-"3Small Eukaryotes"

richEukN2$X.SampleID<-rownames(richEukN2)
richEukN2$X.SampleID<-gsub(pattern = "N", replace = "S", x = richEukN2$X.SampleID)
richEukN3<-merge(richEukN2,biogeo6,"X.SampleID")
richEukN3$type<-"4Soil Microfauna"


richdata<-rbind(richBac3,richITS3,richEukS3,richEukN3)

richmeans<-richdata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(PD),se=std.error(PD))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/diversitybysuccessionalstagechao1.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(richmeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
#dev.off()

# anova(mb<-lm(PD~lomehi,data=richBac3))
# anova(mi<-lm(PD~lomehi,data=richITS3))
# anova(ms<-lm(PD~lomehi,data=richEukS3))
# anova(mn<-lm(PD~lomehi,data=richEukN3))
# 
# summary(mb <- aov(PD~lomehi, data = richBac3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# plot(TukeyHSD(mb, "lomehi"))
# summary(mb <- aov(PD~lomehi, data = richITS3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# summary(mb <- aov(PD~lomehi, data = richEukS3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# summary(mb <- aov(PD~lomehi, data = richEukN3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)

richdatadist<-cbind(richdata$X,richdata$Y,richdata$elevation)
hmiscDISTm<-as.matrix(dist(richdatadist))

##### models with lomehi categorical and spatial autocorrelation #####
#adding spatial autocorrelation, I will not use a nugget b/c it is not significant and b/c the distribution modeling did not use a nugget
richBac3$lomehi<-as.factor(richBac3$lomehi)
m0<-gls(PD~lomehi,data=richBac3)
m1<-gls(PD~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=richBac3)
#m2<-gls(PD~lomehi,correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richBac3)
# anova(m0,m1)
# anova(m0)
anova(m1)
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)}
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

richITS3$lomehi<-as.factor(richITS3$lomehi)
m0<-gls(PD~lomehi,data=richITS3)
m1<-gls(PD~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=richITS3)
# anova(m0,m1)
# anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

richEukS3$lomehi<-as.factor(richEukS3$lomehi)
m0<-gls(PD~lomehi,data=richEukS3)
m1<-gls(PD~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=richEukS3)
# anova(m0,m1)
# anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

richEukN3$lomehi<-as.factor(richEukN3$lomehi)
m0<-gls(PD~lomehi,data=richEukN3)
m1<-gls(PD~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=richEukN3)
# anova(m0,m1)
# anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))





###### Succession as continuous with autocorrelation #####
#I need to log transform, otherwise it is very skewed
ggplot(richdata,aes(x=log(pca1+1),y=PD))+# as.numeric(fert),color=species
  labs(x="Succession",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")

ggplot(richdata,aes(x=pca1,y=PD))+# as.numeric(fert),color=species
  labs(x="Succession",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")

m0<-gls(PD~log(pca1+1),data=richBac3)
m1<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richBac3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richBac3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(PD~log(pca1+1),data=richITS3)
m1<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richITS3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richITS3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(PD~log(pca1+1),data=richEukS3)
m1<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richEukS3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richEukS3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(PD~log(pca1+1),data=richEukN3)
m1<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richEukN3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richEukN3)
anova(m0,m1)
anova(m0)
anova(m1)





##### Species richness #####

richBac2$X.SampleID<-rownames(richBac2)
richBac3<-merge(richBac2,biogeo6,"X.SampleID")
richBac3$type<-"1Bacteria"

richITS2$X.SampleID<-rownames(richITS2)
richITS3<-merge(richITS2,biogeo6,"X.SampleID")
richITS3$type<-"2Fungi"

richEukS2$X.SampleID<-rownames(richEukS2)
richEukS3<-merge(richEukS2,biogeo6,"X.SampleID")
richEukS3$type<-"3Small Eukaryotes"

richEukN2$X.SampleID<-rownames(richEukN2)
richEukN2$X.SampleID<-gsub(pattern = "N", replace = "S", x = richEukN2$X.SampleID)
richEukN3<-merge(richEukN2,biogeo6,"X.SampleID")
richEukN3$type<-"4Soil Microfauna"

richdata<-rbind(richBac3,richITS3,richEukS3,richEukN3)

richmeans<-richdata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(Chao1),se=std.error(Chao1))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/richnessbysuccessionalstagechao1.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(richmeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Taxonomic Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()

# anova(lm(Chao1~lomehi,data=richBac3))
# anova(lm(Chao1~lomehi,data=richITS3))
# anova(lm(Chao1~lomehi,data=richEukS3))
# anova(lm(Chao1~lomehi,data=richEukN3))
# 
# summary(mb <- aov(Chao1~lomehi, data = richBac3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# plot(TukeyHSD(mb, "lomehi"))
# summary(mb <- aov(Chao1~lomehi, data = richITS3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# summary(mb <- aov(Chao1~lomehi, data = richEukS3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# summary(mb <- aov(Chao1~lomehi, data = richEukN3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)


ggplot(richdata,aes(x=log(pca1+1),y=Chao1))+# as.numeric(fert),color=species
  labs(x="Succession",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")

m0<-gls(Chao1~log(pca1+1),data=richBac3)
m1<-gls(Chao1~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richBac3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(Chao1~log(pca1+1),data=richITS3)
m1<-gls(Chao1~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richITS3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richITS3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(Chao1~log(pca1+1),data=richEukS3)
m1<-gls(Chao1~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richEukS3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richEukS3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(Chao1~log(pca1+1),data=richEukN3)
m1<-gls(Chao1~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richEukN3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richEukN3)
anova(m0,m1)
anova(m0)
anova(m1)





##### Evenness #####

richBac2$Evenness<-vegan::diversity(datBacS5cotu[,-c(1:33)])/log(specnumber(datBacS5cotu[,-c(1:33)]))
richBac3<-merge(richBac2,biogeo6,"X.SampleID")
richBac3$type<-"1Bacteria"

richITS2$Evenness<-vegan::diversity(datITSS5cotu[,-c(1:33)])/log(specnumber(datITSS5cotu[,-c(1:33)]))
richITS3<-merge(richITS2,biogeo6,"X.SampleID")
richITS3$type<-"2Fungi"

richEukS2$Evenness<-vegan::diversity(datEukS5cotu[,-c(1:33)])/log(specnumber(datEukS5cotu[,-c(1:33)]))
richEukS3<-merge(richEukS2,biogeo6,"X.SampleID")
richEukS3$type<-"3Small Eukaryotes"

richEukN2$Evenness<-vegan::diversity(datEukN5cotu[,-c(1:33)])/log(specnumber(datEukN5cotu[,-c(1:33)]))
richEukN3<-merge(richEukN2,biogeo6,"X.SampleID")
richEukN3$type<-"4Soil Mesofauna"

richdata<-rbind(richBac3,richITS3,richEukS3,richEukN3)

richmeans<-richdata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(Evenness,na.rm=T),se=std.error(Evenness,na.rm=T))
richmeans$lomehi<-factor(richmeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforMolEcolSubmission/evennessbysuccessionalstage.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(richmeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Evenness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
# dev.off()


# anova(lm(Evenness~lomehi,data=richBac3))
# anova(lm(Evenness~lomehi,data=richITS3))
# anova(lm(Evenness~lomehi,data=richEukS3))
# anova(lm(Evenness~lomehi,data=richEukN3))
# 
# summary(mb <- aov(Evenness~lomehi, data = richBac3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# plot(TukeyHSD(mb, "lomehi"))
# summary(mb <- aov(Evenness~lomehi, data = richITS3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# summary(mb <- aov(Evenness~lomehi, data = richEukS3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)
# summary(mb <- aov(Evenness~lomehi, data = richEukN3))
# TukeyHSD(mb, "lomehi", ordered = TRUE)

richBac3$lomehi<-as.factor(richBac3$lomehi)
m0<-gls(Evenness~lomehi,data=richBac3)
m1<-gls(Evenness~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=richBac3)
#m2<-gls(PD~lomehi,correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richBac3)
# anova(m0,m1)
# anova(m0)
anova(m1)
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)}
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

richITS3$lomehi<-as.factor(richITS3$lomehi)
m0<-gls(Evenness~lomehi,data=richITS3)
m1<-gls(Evenness~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=richITS3)
# anova(m0,m1)
# anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

richEukS3$lomehi<-as.factor(richEukS3$lomehi)
m0<-gls(Evenness~lomehi,data=richEukS3)
m1<-gls(Evenness~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=richEukS3)
# anova(m0,m1)
# anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

richEukN3$lomehi<-as.factor(richEukN3$lomehi)
m0<-gls(Evenness~lomehi,data=richEukN3,na.action = na.omit)
m1<-gls(Evenness~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=richEukN3,na.action = na.omit)
# anova(m0,m1)
# anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))


ggplot(richdata,aes(x=log(pca1+1),y=Evenness))+# as.numeric(fert),color=species
  labs(x="Succession",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")

m0<-gls(Evenness~log(pca1+1),data=richBac3)
m1<-gls(Evenness~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richBac3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(Evenness~log(pca1+1),data=richITS3)
m1<-gls(Evenness~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richITS3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richITS3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(Evenness~log(pca1+1),data=richEukS3)
m1<-gls(Evenness~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richEukS3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richEukS3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(Evenness~log(pca1+1),data=richEukN3,na.action = na.omit)
m1<-gls(Evenness~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richEukN3,na.action = na.omit)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richEukN3)
anova(m0,m1)
anova(m0)
anova(m1)



###### Rarity ######
#Rarity was calculated as the proportion of species that had relative abundances less than 1/S (Camargo 1992). this doesn't really work b/c it standardizes over spcies richness, so for rich plots the cutoff is lower than for non rich plots so they end up having the same rarity which makes no sense. so I should choose a consistent relativ abundance and stick with that across all plots.

richEukS2
richEukN2
richBac2
richITS2

datEukS5cotu
datEukN5cotu
datBacS5cotu
datITSS5cotu


#bac
richBac2$X.SampleID<-rownames(richBac2)
richBac3<-merge(richBac2,biogeo6,"X.SampleID")
datBacS5otu2<-merge(biogeo6,datBacS5otu,"X.SampleID")
row.names(datBacS5otu2)<-datBacS5otu2$X.SampleID
datBacSp<-datBacS5otu2[,87:dim(datBacS5otu2)[2]]
ind<-which(colSums(datBacSp)>0)
datBacSp2<-datBacSp[,ind]

rarityBac<-data.frame("X.SampleID"=row.names(datBacSp2),rarity=NA)
for(i in 1:75){
  tempr<-1/mean(richBac3$SR)
  #tempr<-.001
  lengthr<-length(which(datBacSp2[i,]>0&datBacSp2[i,]<tempr))
  rarityBac$rarity[i]<-lengthr/richBac3$SR[i]
}
rarityBac2<-merge(rarityBac,biogeo6,"X.SampleID")
rarityBac2$type<-"1Bacteria"

rarityBac2%>%
  mutate(lomehi=factor(lomehi,levels=c('lo','me','hi')))%>%
  group_by(lomehi)%>%
  summarise(mean=mean(rarity),se=std.error(rarity))


#ITS
richITS2$X.SampleID<-rownames(richITS2)
richITS3<-merge(richITS2,biogeo6,"X.SampleID")
datITSS5otu2<-merge(biogeo6,datITSS5otu,"X.SampleID")
row.names(datITSS5otu2)<-datITSS5otu2$X.SampleID
datITSSp<-datITSS5otu2[,87:dim(datITSS5otu2)[2]]
ind<-which(colSums(datITSSp)>0)
datITSSp2<-datITSSp[,ind]

rarityITS<-data.frame("X.SampleID"=row.names(datITSSp2),rarity=NA)
for(i in 1:75){
  tempr<-1/mean(richITS3$SR)
  lengthr<-length(which(datITSSp2[i,]>0&datITSSp2[i,]<tempr))
  rarityITS$rarity[i]<-lengthr/richITS3$SR[i]
}
rarityITS2<-merge(rarityITS,biogeo6,"X.SampleID")
rarityITS2$type<-"2Fungi"

rarityITS2%>%
  mutate(lomehi=factor(lomehi,levels=c('lo','me','hi')))%>%
  group_by(lomehi)%>%
  summarise(mean=mean(rarity),se=std.error(rarity))

#euk S
richEukS2$X.SampleID<-rownames(richEukS2)
richEukS3<-merge(richEukS2,biogeo6,"X.SampleID")
datEukS5otu2<-merge(biogeo6,datEukS5otu,"X.SampleID")
row.names(datEukS5otu2)<-datEukS5otu2$X.SampleID
datEukSp<-datEukS5otu2[,87:dim(datEukS5otu2)[2]]
ind<-which(colSums(datEukSp)>0)
datEukSp2<-datEukSp[,ind]

rarityEuk<-data.frame("X.SampleID"=row.names(datEukSp2),rarity=NA)
for(i in 1:75){
  tempr<-1/mean(richEukS3$SR)
  #tempr<-.01
  lengthr<-length(which(datEukSp2[i,]>0&datEukSp2[i,]<tempr))
  rarityEuk$rarity[i]<-lengthr/richEukS3$SR[i]
}
rarityEuk2<-merge(rarityEuk,biogeo6,"X.SampleID")
rarityEuk2$type<-"3Small Eukaryotes"

rarityEuk2%>%
  mutate(lomehi=factor(lomehi,levels=c('lo','me','hi')))%>%
  group_by(lomehi)%>%
  summarise(mean=mean(rarity),se=std.error(rarity))

#euk N
richEukN2$X.SampleID<-rownames(richEukN2)
richEukN2$X.SampleID<-gsub(pattern = "N", replace = "S", x = richEukN2$X.SampleID)
richEukN3<-merge(richEukN2,biogeo6,"X.SampleID")
datEukN6otu<-datEukN5otu
datEukN6otu$X.SampleID<-gsub(pattern = "N", replace = "S", x = datEukN6otu$X.SampleID)
datEukN5otu2<-merge(biogeo6,datEukN6otu,"X.SampleID")
row.names(datEukN5otu2)<-datEukN5otu2$X.SampleID
datEukNSp<-datEukN5otu2[,87:dim(datEukN5otu2)[2]]
ind<-which(colSums(datEukNSp)>0)
datEukNSp2<-datEukNSp[,ind]

rarityEukN<-data.frame("X.SampleID"=row.names(datEukNSp2),rarity=NA)
for(i in 1:75){
  tempr<-1/mean(richEukN3$SR)
  #tempr<-.01
  lengthr<-length(which(datEukNSp2[i,]>0&datEukNSp2[i,]<tempr))
  rarityEukN$rarity[i]<-lengthr/richEukN3$SR[i]
}
rarityEukN2<-merge(rarityEukN,biogeo6,"X.SampleID")
rarityEukN2$type<-"4Soil Microfauna"

rarityEukN2%>%
  mutate(lomehi=factor(lomehi,levels=c('lo','me','hi')))%>%
  group_by(lomehi)%>%
  summarise(mean=mean(rarity),se=std.error(rarity))


raritydata<-rbind(rarityBac2,rarityITS2,rarityEuk2,rarityEukN2)

raritymeans<-raritydata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(rarity),se=std.error(rarity))
raritymeans$lomehi<-factor(raritymeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforFrontiersSubmission/raritybysuccessionalstage.pdf",width=3.386,height=3.386) 
ggplot(raritymeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Proportion rare taxa")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
#dev.off()

rarityBac2$lomehi<-factor(rarityBac2$lomehi)
m0<-gls(rarity~lomehi,data=rarityBac2)
m1<-gls(rarity~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=rarityBac2)
# anova(m0,m1)
# anova(m0)
anova(m1)
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)}
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

rarityITS2$lomehi<-factor(rarityITS2$lomehi)
m0<-gls(rarity~lomehi,data=rarityITS2)
m1<-gls(rarity~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=rarityITS2)
anova(m0,m1)
anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

rarityEuk2$lomehi<-factor(rarityEuk2$lomehi)
m0<-gls(rarity~lomehi,data=rarityEuk2)
m1<-gls(rarity~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=rarityEuk2)
anova(m0,m1)
anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))

rarityEukN2$lomehi<-factor(rarityEukN2$lomehi)
m0<-gls(rarity~lomehi,data=rarityEukN2)
m1<-gls(rarity~lomehi,correlation=corSpher(form = ~ X + Y + elevation),data=rarityEukN2)
anova(m0,m1)
anova(m0)
anova(m1)
summary(glht(m1, linfct = mcp(lomehi = "Tukey")))



m0<-gls(rarity~log(pca1+1),data=rarityBac2)
m1<-gls(rarity~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richBac3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(rarity~log(pca1+1),data=richITS3)
m1<-gls(rarity~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richITS3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richITS3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(rarity~log(pca1+1),data=richEukS3)
m1<-gls(rarity~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richEukS3)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richEukS3)
anova(m0,m1)
anova(m0)
anova(m1)

m0<-gls(rarity~log(pca1+1),data=richEukN3,na.action = na.omit)
m1<-gls(rarity~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation),data=richEukN3,na.action = na.omit)
#m2<-gls(PD~log(pca1+1),correlation=corSpher(form = ~ X + Y + elevation,nugget=T),data=richEukN3)
anova(m0,m1)
anova(m0)
anova(m1)






###### Frequency #####

#bac
# richBac2$X.SampleID<-rownames(richBac2)
# richBac3<-merge(richBac2,biogeo6,"X.SampleID")
datBacS5otu2<-merge(biogeo6,datBacS5otu,"X.SampleID")
row.names(datBacS5otu2)<-datBacS5otu2$X.SampleID
datBacSp<-datBacS5otu2[,87:dim(datBacS5otu2)[2]]
ind<-which(colSums(datBacSp)>0)
datBacSp2<-datBacSp[,ind]

ind<-which(datBacS5otu2$lomehi.x=="lo")
datBacSp2lo<-datBacSp2[ind,]
ind<-which(colSums(datBacSp2lo)>0)
datBacSp3lo<-datBacSp2lo[,ind]
mean(colSums(datBacSp3lo>0))
#percent taxa that we are modeling in network analysis
length(which(colSums(datBacSp3lo>0)>11))/dim(datBacSp3lo)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datBacSp3lo>0)<3))/dim(datBacSp3lo)[2]

ind<-which(datBacS5otu2$lomehi.x=="me")
datBacSp2me<-datBacSp2[ind,]
ind<-which(colSums(datBacSp2me)>0)
datBacSp3me<-datBacSp2me[,ind]
mean(colSums(datBacSp3me>0))
length(which(colSums(datBacSp3me>0)>11))/dim(datBacSp3me)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datBacSp3me>0)<3))/dim(datBacSp3me)[2]

ind<-which(datBacS5otu2$lomehi.x=="hi")
datBacSp2hi<-datBacSp2[ind,]
ind<-which(colSums(datBacSp2hi)>0)
datBacSp3hi<-datBacSp2hi[,ind]
mean(colSums(datBacSp3hi>0))
length(which(colSums(datBacSp3hi>0)>11))/dim(datBacSp3hi)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datBacSp3hi>0)<3))/dim(datBacSp3hi)[2]

freqdatabaclo<-data.frame(freq=colSums(datBacSp3lo>0),lomehi="lo")
freqdatabacme<-data.frame(freq=colSums(datBacSp3me>0),lomehi="me")
freqdatabachi<-data.frame(freq=colSums(datBacSp3hi>0),lomehi="hi")
freqdatabac<-rbind(freqdatabaclo,freqdatabacme,freqdatabachi)
freqdatabac$type<-"1Bacteria"


#its
datITSS5otu2<-merge(biogeo6,datITSS5otu,"X.SampleID")
row.names(datITSS5otu2)<-datITSS5otu2$X.SampleID
datITSSp<-datITSS5otu2[,87:dim(datITSS5otu2)[2]]
ind<-which(colSums(datITSSp)>0)
datITSSp2<-datITSSp[,ind]

ind<-which(datITSS5otu2$lomehi.x=="lo")
datITSSp2lo<-datITSSp2[ind,]
ind<-which(colSums(datITSSp2lo)>0)
datITSSp3lo<-datITSSp2lo[,ind]
mean(colSums(datITSSp3lo>0))
#percent taxa that we are modeling in network analysis
length(which(colSums(datITSSp3lo>0)>11))/dim(datITSSp3lo)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datITSSp3lo>0)<3))/dim(datITSSp3lo)[2]

ind<-which(datITSS5otu2$lomehi.x=="me")
datITSSp2me<-datITSSp2[ind,]
ind<-which(colSums(datITSSp2me)>0)
datITSSp3me<-datITSSp2me[,ind]
mean(colSums(datITSSp3me>0))
length(which(colSums(datITSSp3me>0)>11))/dim(datITSSp3me)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datITSSp3me>0)<3))/dim(datITSSp3me)[2]

ind<-which(datITSS5otu2$lomehi.x=="hi")
datITSSp2hi<-datITSSp2[ind,]
ind<-which(colSums(datITSSp2hi)>0)
datITSSp3hi<-datITSSp2hi[,ind]
mean(colSums(datITSSp3hi>0))
length(which(colSums(datITSSp3hi>0)>11))/dim(datITSSp3hi)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datITSSp3hi>0)<3))/dim(datITSSp3hi)[2]

freqdataITSlo<-data.frame(freq=colSums(datITSSp3lo>0),lomehi="lo")
freqdataITSme<-data.frame(freq=colSums(datITSSp3me>0),lomehi="me")
freqdataITShi<-data.frame(freq=colSums(datITSSp3hi>0),lomehi="hi")
freqdataITS<-rbind(freqdataITSlo,freqdataITSme,freqdataITShi)
freqdataITS$type<-"2Fungi"

#eukS
datEukS5otu2<-merge(biogeo6,datEukS5otu,"X.SampleID")
row.names(datEukS5otu2)<-datEukS5otu2$X.SampleID
datEukSp<-datEukS5otu2[,87:dim(datEukS5otu2)[2]]
ind<-which(colSums(datEukSp)>0)
datEukSp2<-datEukSp[,ind]

ind<-which(datEukS5otu2$lomehi.x=="lo")
datEukSp2lo<-datEukSp2[ind,]
ind<-which(colSums(datEukSp2lo)>0)
datEukSp3lo<-datEukSp2lo[,ind]
mean(colSums(datEukSp3lo>0))
#percent taxa that we are modeling in network analysis
length(which(colSums(datEukSp3lo>0)>11))/dim(datEukSp3lo)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datEukSp3lo>0)<3))/dim(datEukSp3lo)[2]

ind<-which(datEukS5otu2$lomehi.x=="me")
datEukSp2me<-datEukSp2[ind,]
ind<-which(colSums(datEukSp2me)>0)
datEukSp3me<-datEukSp2me[,ind]
mean(colSums(datEukSp3me>0))
length(which(colSums(datEukSp3me>0)>11))/dim(datEukSp3me)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datEukSp3me>0)<3))/dim(datEukSp3me)[2]

ind<-which(datEukS5otu2$lomehi.x=="hi")
datEukSp2hi<-datEukSp2[ind,]
ind<-which(colSums(datEukSp2hi)>0)
datEukSp3hi<-datEukSp2hi[,ind]
mean(colSums(datEukSp3hi>0))
length(which(colSums(datEukSp3hi>0)>11))/dim(datEukSp3hi)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datEukSp3hi>0)<3))/dim(datEukSp3hi)[2]

freqdataEuklo<-data.frame(freq=colSums(datEukSp3lo>0),lomehi="lo")
freqdataEukme<-data.frame(freq=colSums(datEukSp3me>0),lomehi="me")
freqdataEukhi<-data.frame(freq=colSums(datEukSp3hi>0),lomehi="hi")
freqdataEuk<-rbind(freqdataEuklo,freqdataEukme,freqdataEukhi)
freqdataEuk$type<-"3Small Eukaryotes"

#eukN
datEukN5otub<-datEukN5otu
datEukN5otub$X.SampleID<-gsub(pattern = "N", replace = "S", x = datEukN5otub$X.SampleID)
datEukN5otu2<-merge(biogeo6,datEukN5otub,"X.SampleID")
row.names(datEukN5otu2)<-datEukN5otu2$X.SampleID
datEukSp<-datEukN5otu2[,87:dim(datEukN5otu2)[2]]
ind<-which(colSums(datEukNSp)>0)
datEukNSp2<-datEukNSp[,ind]

ind<-which(datEukN5otu2$lomehi.x=="lo")
datEukNSp2lo<-datEukNSp2[ind,]
ind<-which(colSums(datEukNSp2lo)>0)
datEukNSp3lo<-datEukNSp2lo[,ind]
mean(colSums(datEukNSp3lo>0))
#percent taxa that we are modeling in network analysis
length(which(colSums(datEukNSp3lo>0)>11))/dim(datEukNSp3lo)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datEukNSp3lo>0)<3))/dim(datEukNSp3lo)[2]

ind<-which(datEukN5otu2$lomehi.x=="me")
datEukNSp2me<-datEukNSp2[ind,]
ind<-which(colSums(datEukNSp2me)>0)
datEukNSp3me<-datEukNSp2me[,ind]
mean(colSums(datEukNSp3me>0))
#percent taxa that we are modeling in network analysis
length(which(colSums(datEukNSp3me>0)>11))/dim(datEukNSp3me)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datEukNSp3me>0)<3))/dim(datEukNSp3me)[2]

ind<-which(datEukN5otu2$lomehi.x=="hi")
datEukNSp2hi<-datEukNSp2[ind,]
ind<-which(colSums(datEukNSp2hi)>0)
datEukNSp3hi<-datEukNSp2hi[,ind]
mean(colSums(datEukNSp3hi>0))
#percent taxa that we are modeling in network analysis
length(which(colSums(datEukNSp3hi>0)>11))/dim(datEukNSp3hi)[2]
#percent taxa that are found in only 1 or 2 plots
length(which(colSums(datEukNSp3hi>0)<3))/dim(datEukNSp3hi)[2]

freqdataEukNlo<-data.frame(freq=colSums(datEukNSp3lo>0),lomehi="lo")
freqdataEukNme<-data.frame(freq=colSums(datEukNSp3me>0),lomehi="me")
freqdataEukNhi<-data.frame(freq=colSums(datEukNSp3hi>0),lomehi="hi")
freqdataEukN<-rbind(freqdataEukNlo,freqdataEukNme,freqdataEukNhi)
freqdataEukN$type<-"4Soil Microfauna"

freqdata<-rbind(freqdatabac,freqdataITS,freqdataEuk,freqdataEukN)

freqmeans<-freqdata%>%
  group_by(type,lomehi)%>%
  summarise(mean=mean(freq,na.rm=T),se=std.error(freq,na.rm=T))
freqmeans$lomehi<-factor(freqmeans$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/Figs/FigsforFrontiersSubmission/freqbysuccessionalstage.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(freqmeans,aes(x=lomehi,y=mean,group=type))+
  labs(x = "",y="Frequency (number of plots)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)+#.5
  #scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
# dev.off()


freqdatabac$taxon<-rownames(freqdatabac)
m0<-lme(freq~lomehi,random=~1|taxon,data=freqdatabac)
m1<-gls(freq~lomehi,data=freqdatabac)
m2<-lm(freq~lomehi,data=freqdatabac)
m3<-glm(freq~lomehi,family="poisson",data=freqdatabac)
m4<-glm(freq~lomehi,family="quasipoisson",data=freqdatabac)
AIC(m2,m3,m4)
drop1(m4,test="F")
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)}
summary(glht(m4, linfct = mcp(lomehi = "Tukey")))

freqdataITS$taxon<-rownames(freqdataITS)
m0<-lme(freq~lomehi,random=~1|taxon,data=freqdataITS)
m1<-gls(freq~lomehi,data=freqdataITS)
anova(m0,m1)
m2<-lm(freq~lomehi,data=freqdataITS)
m3<-glm(freq~lomehi,family="poisson",data=freqdataITS)
m4<-glm(freq~lomehi,family="quasipoisson",data=freqdataITS)
AIC(m2,m3,m4)
drop1(m4,test="F")
summary(glht(m4, linfct = mcp(lomehi = "Tukey")))

freqdataEuk$taxon<-rownames(freqdataEuk)
m0<-lme(freq~lomehi,random=~1|taxon,data=freqdataEuk)
m1<-gls(freq~lomehi,data=freqdataEuk)
anova(m0,m1)
m2<-lm(freq~lomehi,data=freqdataEuk)
m3<-glm(freq~lomehi,family="poisson",data=freqdataEuk)
m4<-glm(freq~lomehi,family="quasipoisson",data=freqdataEuk)
AIC(m2,m3,m4)
drop1(m4,test="F")
summary(glht(m4, linfct = mcp(lomehi = "Tukey")))

freqdataEukN$taxon<-rownames(freqdataEukN)
m0<-lme(freq~lomehi,random=~1|taxon,data=freqdataEukN)
m1<-gls(freq~lomehi,data=freqdataEukN)
anova(m0,m1)
m2<-lm(freq~lomehi,data=freqdataEukN)
m3<-glm(freq~lomehi,family="poisson",data=freqdataEukN)
m4<-glm(freq~lomehi,family="quasipoisson",data=freqdataEukN)
AIC(m2,m3,m4)
drop1(m4,test="F")
summary(glht(m4, linfct = mcp(lomehi = "Tukey")))













###### Old code when testing DADA2 with different parameters (none of them made a difference) ######

biogeo8<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/biogeo8.csv")

mapBac<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/515BC_Niwot_20072015_All_MapFilenewlomehinoN472015.txt")
mapITS<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITSsingle/ITS_Niwot_20072015_All_MapFilenewlomehi.txt")

ps <- phyloseq(otu_table(seqtab.nochimb, taxa_are_rows=FALSE), 
               sample_data(mapBac))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(mapITS))

ps2<-ps%>%
  subset_samples(SampleType=="soil"&year==2015)%>%
  subset_samples(Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)%>%
  filter_taxa(function(x) sum(x) > (1), prune=T)%>%##even if there aren't singletons, it prunes them which is nice 
  rarefy_even_depth(sample.size=1023,rngseed=10,replace=F) #%>% for bac 6919 for its 1056
#transform_sample_counts(function(x) x/sum(x) )
sort(sample_sums(ps2))

plot_richness(ps2, x="VascPlant_Div", measures=c("Chao1", "Observed"))

ps3<-cbind(sample_data(ps2),otu_table(ps2))
ps3$Sample_name<-as.numeric(as.character(ps3$Sample_name))

ps3$richness<-rowSums(ps3[,32:dim(ps3)[2]]>0)
m1<-aggregate.data.frame(ps3$richness,by=list(ps3$lomehi),mean)
se<-aggregate.data.frame(ps3$richness,by=list(ps3$lomehi),std.error)
plot(log(ps3$VascPlant_Dens+1),ps3$richness)
abline(lm(ps3$richness~log(ps3$VascPlant_Dens+1)))
summary(lm(ps3$richness~log(ps3$VascPlant_Dens+1)))
plot(ps3$VascPlant_Div,ps3$richness)
abline(lm(ps3$richness~ps3$VascPlant_Div))

ps3means<-ps3%>%
  group_by(lomehi)%>%
  summarise(mean=mean(richness),se=std.error(richness))
ps3means$lomehi<-factor(ps3means$lomehi,levels=c("lo","me","hi"))

#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/diversitybysuccessionalstage.pdf",width=3.386,height=3.386) #width=3.386 or 7
ggplot(ps3means,aes(x=lomehi,y=mean))+
  labs(x = "",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.4)+#.5
  geom_point(size=1.5)+#2
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.4)






