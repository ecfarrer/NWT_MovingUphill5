#Submitting to Genbank

#go here to get soil and environment codes
http://bioportal.bioontology.org/ontologies/ENVO/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_00005741
http://bioportal.bioontology.org/ontologies/ENVO/?p=classes&conceptid=http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FENVO_00000112


#to download and install aspera go https://www.ncbi.nlm.nih.gov/sra/docs/submitfiles/ and then to https://downloads.asperasoft.com/en/downloads/8?list
#do the install as usual.
#when you are on the page to upload the files a popup happens that asks you if you want to initiate the aspera plug in and I said yes. then you just do everything by command line. I can't remember if you actually have to install the plugin in chrome. I don't think I did it for my laptop, but I did for my desktop b/c I initially went to a different window in the aspera page to try to download it and a window popped up saying I didn't have the chrome plugin.


#####Bacteria#### 

#the demultiplexed files were in /Bact/exported_demux (which is the same as what is in /Bact/exported_demuxdada2, they were copied over to the exported_demuxdada2 folder to complete the dada2 bioinformatics, which actually was not used in the final analysis. instead I used the qiime2 pipeline)

#I created this directly from the tsv file used in QIIME2
genbank16S<-read.csv("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/515BC_Niwot_20072015_All_MapFilenewlomehinoN472015FORGENBANKSUBMISSION.csv")

head(genbank16S)

#Coordinates
snowdepth<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/King_SnowDepth.csv")
head(snowdepth)
colnames(snowdepth)[1]<-"Sample_name"


library(rgdal)
latlong16S <- SpatialPoints(snowdepth[,c("EASTING","NORTHING")], proj4string=CRS("+proj=utm +zone=13 +datum=WGS84"))
latlong16S2 <- spTransform(latlong16S, CRS("+proj=longlat +datum=WGS84"))
snowdepth<-cbind(snowdepth,latlong16S2@coords)
colnames(snowdepth)[45:46]<-c("lon","lat")
head(snowdepth)
snowdepth2<-snowdepth[,c(1,45:46)]
snowdepth2$lon<-round(snowdepth2$lon,5)
snowdepth2$lat<-round(snowdepth2$lat,5)

genbank16S2<-merge(genbank16S,snowdepth2)
head(genbank16S2)

#format: 38.98 N 77.11 W

genbank16S2$latlon<-paste(genbank16S2$lat,"N",-genbank16S2$lon,"W",sep=" ")

genbank16S2$latlon
genbank16S2$elevation

output16S<-cbind(X.SampleID=as.character(genbank16S2$X.SampleID),latlon=genbank16S2$latlon,elevation=genbank16S2$elevation,date=as.character(genbank16S2$collection_date))
  
write.csv(output16S,"~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/output16S.csv")             



#To rename output demux files in unix for bacteria
#to remove the _001
for f in S*; do mv -n "$f" "${f/_001/}";done;ls -l

#add 16S
for f in S*; do mv -n "$f" "${f/_L001/MII16S}";done;ls -l

#replace second _ with period
for f in *; do mv "$f" "$(echo "$f" | sed 's/_/./2')"; done;ls -l

# delete everything between _ and MII with period
for f in *; do mv "$f" "$(echo "$f" | sed 's/_.*MII/./')"; done;ls -l


#for ftp uploads
#first navigate to your folder
#I needed to intsall the ftp function: brew install inetutils
ftp
open ftp-private.ncbi.nlm.nih.gov
subftp
w4pYB9VQ
cd uploads/efarrer_tulane.edu_zJzARsfJ
mkdir 16S_2015
cd 16S_2015
mput *.fastq.gz
#you need to push enter after each entry
 
#the ftp upload did not work, said files were corrupted.
#trying aspera
#this kept failing, so I tried lots of different things. on thing I did was add 130.14.29.0/24 and 130.14.250.0/24 as trustd hosts in aspera connect. I have no idea if this is necessary.
/Users/farrer/Applications/Aspera\ Connect.app/Contents/Resources/ascp -i /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/aspera.openssh -QT -l100m -k1 -d /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Bact/GenBankSRASubmission16S subasp@upload.ncbi.nlm.nih.gov:uploads/efarrer_tulane.edu_chwkHMsu




##### ITS #####
#the demultiplexed files were in /ITS/exported_demux (which is the same as what is in /ITS/exported_demuxdada2, they were copied over to the exported_demuxdada2 folder to complete the dada2 bioinformatics)

#To rename output demux files in unix for ITS
#to remove the _001
for f in S*; do mv -n "$f" "${f/_001/}";done;ls -l

#add ITS
for f in S*; do mv -n "$f" "${f/_L001/MOOITS}";done;ls -l

#replace second _ with period
for f in *; do mv "$f" "$(echo "$f" | sed 's/_/./2')"; done;ls -l

# delete everything between _ and MII with period
for f in *; do mv "$f" "$(echo "$f" | sed 's/_.*MOO/./')"; done;ls -l

#Look inside a gz file
zcat < S.159.2015.ITS.R1.fastq.gz | tail

#commandline aspera
/Users/farrer/Applications/Aspera\ Connect.app/Contents/Resources/ascp -i /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/aspera.openssh -QT -l100m -k1 -d /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/ITS/GenBankSRASubmissionITS subasp@upload.ncbi.nlm.nih.gov:uploads/efarrer_tulane.edu_chwkHMsu






##### 18S #####
#the demultiplexed files were in /Euks/exported_demux. there are no other demultiplexed files b/c I did paired end analysis for euks from the start 

#To rename output demux files in unix for ITS
#to remove the _001
for f in *; do mv -n "$f" "${f/_001/}";done;ls -l

#add 18S
for f in *; do mv -n "$f" "${f/_L001/MOO18S}";done;ls -l

#replace second _ with period
for f in *; do mv "$f" "$(echo "$f" | sed 's/_/./2')"; done;ls -l

# delete everything between _ and MII with period
for f in *; do mv "$f" "$(echo "$f" | sed 's/_.*MOO/./')"; done;ls -l

#commandline aspera
/Users/farrer/Applications/Aspera\ Connect.app/Contents/Resources/ascp -i /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Euks/aspera.openssh -QT -l100m -k1 -d /Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/QIIME2/Euks/GenBankSRASubmission18S subasp@upload.ncbi.nlm.nih.gov:uploads/efarrer_tulane.edu_chwkHMsu


