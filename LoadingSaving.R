#Order of running files:
#DataCleaning.R: to clean microbe and plant data
#SuccessionalStagePCA.R: to incorporate biogeochemistry and do a PCA for successional stage
#boral: doing a little more data reorganizing, calculating clr, running hirarchical models, graphing networks
#Randomizations.R: 


#Loading/saving/packages needed

#if vetor memory exhausted happens do this in terminal
cd ~
touch .Renviron #create this file if it doesn't already exist
open .Renviron
#then write this in the file that opens, need to play with this to see what crashes R. 700 was too much. 70 on my laptop seems to work ok, maybe a little too much but didn't crash, 70 did crash on my laptop, trying 35 on laptop (now trying 40 since on loading envronmnts I got th error), 70 on imac seems good
R_MAX_VSIZE=70Gb 
http://r.789695.n4.nabble.com/R-3-5-0-vector-memory-exhausted-error-on-readBin-td4750237.html
http://btibert3.github.io/2015/12/08/Environment-Variables-in-Rstudio-on-Mac.html


#enviroments:

#These are the same enviroments from the last time, nothing has changed
MovingUphill4_WorkspaceITSBactTrials.Rdata #dada2 trials in R when I was testing truncation and trimming for ITS and bacteria
MovingUphill4_WorkspaceITSbioinformatics.Rdata
MovingUphill4_WorkspaceDataCleaning.Rdata
MovingUphill4_WorkspaceDataCleaningOutput.Rdata #just the 15 output files I need for downstream analysis, all intermediate files deleted from env
MovingUphill4_WorkspaceAnalysisNetworkTrials.Rdata #all the different trials I ran when deciding the parameters for the networks
MovingUphill4_WorkspaceSimulations.Rdata # simulations based on Faust 2015

#These environments are where the new analyses start
MovingUphill5_SubsettingPCArange.Rdata #subsetting low and hi at the same range of pca axis 1
MovingUphill5_WorkspaceSubsetting.Rdata #subsetting hi plots 10 times 273 of the 306 species to make it similar to the number of taxa going in to the network analysis
MovingUphill5_WorkspaceRandomizations.Rdata #randomizing my data to detect false positives
MovingUphill5_Krona.Rdata #code for krona plots
MovingUphill5_WorkspaceTrials1.Rdata #for boral models
MovingUphill5_WorkspaceTrials2.Rdata #for boral models

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/MovingUphill5_WorkspaceTrials1.Rdata")  # 

load("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/FiguresStats/kingdata/MovingUphill5_WorkspaceTrials1.Rdata") 

#rm(list=setdiff(ls(), c("fit.lolv4occ9exp4","rescor.lolv4occ9exp4")))


#for data cleaning

#for installing phyloseq
#biocLite is no longer, now use BiocManager
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

library(phyloseq)
#packageVersion("phyloseq")
library(picante) #for phylogenetic diversity

#### using DADA2 for further processing #####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")

#devtools::install_github("benjjneb/dada2") #to install most recent version when using R3.4, this installs but gives an error, ugh, not dada2 is not available for R3.4, maybe need to do above to get the correct version if it is still available

library(dada2)
packageVersion("dada2") #1.4.0 does not have trimRight
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
# BiocManager::install("decontam")

library(decontam)

#for cooccurrence networks
#library(foreach)
#library(doParallel)

#to install new versions of HMSC
#library(devtools)
#install_github('guiblanchet/HMSC') #takes a long time, 30 min?

#for krona plots
#BiocManager::install("cpauvert/psadd")
library(psadd)

library(HMSC)

library(vegan)
library(corrplot)
library(circlize)
library(Hmisc)
library(boral)
library(Matrix)

#for plotting
library(igraph)
#library(fdrtool)
library(ggplot2)
library(grid) #for unit function in ggplot2 for legend 

library(vegan)
library(multcomp)

#for rarefaction curves
#install.packages("remotes")
#remotes::install_github("gauravsk/ranacapa")
library(ranacapa)

#for network stats
library(NetIndices)

#for manipulating datasets for plotting 
library(tidyr)
library(dplyr)
library(plotrix)

#for doing permutation test, and imputing zeros
library(zCompositions)
#library(combinat)
#library(coin)

#for boral
library(boral) #Need version 0.7 or later, available on CRAN.
library(Matrix)
library(rcompanion) #for nagelkerke R2 on gls models

#for simulations
library(dirmult)
library(HMP)
library(vegan)

#detach(package:igraph)
#sessionInfo()

#extra not needed
#library(reshape)
#library(plotrix)
#library(Kendall)


#library(data.table)
#library(BiodiversityR) #this requires X11 and takes a while to load, you need to close the window that it opens in rcommander


