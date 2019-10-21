# to make a Krona plot for bacteria, eukS and eukN, with lo me hi as the grouping variable

datBacS5
datITSS5
datEukS5
datEukN5

#merge correct lo me hi column from biogeo6 (and reduce plots to only the 75)

biogeo6$lomehi

datBacS6<-datBacS5%>%
  subset_samples(X.SampleID%in%biogeo6$X.SampleID)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

sample_data(datBacS6)$lomehi<-biogeo6$lomehi
sample_data(datBacS6)$SuccessionalStage<-ifelse(sample_data(datBacS6)$lomehi=="lo","Early",ifelse(sample_data(datBacS6)$lomehi=="me","Mid","Late"))

plot_krona(datBacS6,"Bacteria-krona", "SuccessionalStage",trim=F)



datITSS6<-datITSS5%>%
  subset_samples(X.SampleID%in%biogeo6$X.SampleID)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

sample_data(datITSS6)$lomehi<-biogeo6$lomehi
sample_data(datITSS6)$SuccessionalStage<-ifelse(sample_data(datITSS6)$lomehi=="lo","Early",ifelse(sample_data(datITSS6)$lomehi=="me","Mid","Late"))

plot_krona(datITSS6,"Fungi-krona", "SuccessionalStage",trim=F)




datEukS6<-datEukS5%>%
  subset_samples(X.SampleID%in%biogeo6$X.SampleID)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

sample_data(datEukS6)$lomehi<-biogeo6$lomehi
sample_data(datEukS6)$SuccessionalStage<-ifelse(sample_data(datEukS6)$lomehi=="lo","Early",ifelse(sample_data(datEukS6)$lomehi=="me","Mid","Late"))
sample_data(datEukS6)$SuccessionalStage<-factor(sample_data(datEukS6)$SuccessionalStage,levels = c("Early","Mid","Late"))

#something was weird with ranknames (rank 12 and 13 were all NAs so it was messing up so I did it by hand)
#rank_names(datEukS6)<-rank_names(datEukS6)[1:11]
#rank_names(datEukS6)[12:13]<-NULL

#plot_krona(datEukS6,"SmallEukaryotes-krona", "SuccessionalStage",trim=T)

physeq<-datEukS6
variable<-"SuccessionalStage"
output<-"SmallEukaryotes-krona"
trim<-FALSE

function (physeq, output, variable, trim = F) 
{
  if (system(command = "which ktImportText", intern = FALSE, 
             ignore.stdout = TRUE)) {
    stop("KronaTools are not installed. Please see https://github.com/marbl/Krona/wiki/KronaTools.")
  }
  if (is.null(tax_table(physeq))) {
    stop("No taxonomy table available.")
  }
  if (!variable %in% colnames(sample_data(physeq))) {
    stop(paste(variable, "is not a variable in the sample data."))
  }
  if (trim == FALSE) {
    spec.char <- grepl(" |\\(|\\)", as(sample_data(physeq), 
                                       "data.frame")[, variable])
    if (sum(spec.char > 0)) {
      message("The following lines contains spaces or brackets.")
      print(paste(which(spec.char)))
      stop("Use trim=TRUE to convert them automatically or convert manually before re-run")
    }
  }
  df <- psmelt(physeq)
  df <- df[, c("Abundance", variable, rank_names(physeq)[1:11])]
  df[, 2] <- gsub(" |\\(|\\)", "", df[, 2])
  df[, 2] <- as.factor(df[, 2])
  dir.create(output)
  for (lvl in levels(df[, 2])) {
    write.table(unique(df[which(df[, 2] == lvl), -2]), file = paste0(output, 
                                                                     "/", lvl, "taxonomy.txt"), sep = "\t", row.names = F, 
                col.names = F, na = "", quote = F)
  }
  krona_args <- paste(output, "/", levels(df[, 2]), "taxonomy.txt,", 
                      levels(df[, 2]), sep = "", collapse = " ")
  output <- paste(output, ".html", sep = "")
  system(paste("ktImportText", krona_args, "-o", output, sep = " "))
  browseURL(output)
}





datEukN5b<-datEukN5

sample_data(datEukN5b)$X.SampleID<-gsub(pattern = "N", replace = "S", x = sample_data(datEukN5b)$X.SampleID)

sample_data(datEukN5b)
datEukN6<-datEukN5b%>%
  subset_samples(X.SampleID%in%biogeo6$X.SampleID)%>%
  filter_taxa(function(x) sum(x) > (0), prune=T)

sample_data(datEukN6)$lomehi<-biogeo6$lomehi
sample_data(datEukN6)$SuccessionalStage<-ifelse(sample_data(datEukN6)$lomehi=="lo","Early",ifelse(sample_data(datEukN6)$lomehi=="me","Mid","Late"))

#plot_krona(datEukN6,"SoilMicrofauna-krona", "SuccessionalStage",trim=F)

physeq<-datEukN6
variable<-"SuccessionalStage"
output<-"SoilMicrofauna-krona"
trim<-FALSE

function (physeq, output, variable, trim = F) 
{
  if (system(command = "which ktImportText", intern = FALSE, 
             ignore.stdout = TRUE)) {
    stop("KronaTools are not installed. Please see https://github.com/marbl/Krona/wiki/KronaTools.")
  }
  if (is.null(tax_table(physeq))) {
    stop("No taxonomy table available.")
  }
  if (!variable %in% colnames(sample_data(physeq))) {
    stop(paste(variable, "is not a variable in the sample data."))
  }
  if (trim == FALSE) {
    spec.char <- grepl(" |\\(|\\)", as(sample_data(physeq), 
                                       "data.frame")[, variable])
    if (sum(spec.char > 0)) {
      message("The following lines contains spaces or brackets.")
      print(paste(which(spec.char)))
      stop("Use trim=TRUE to convert them automatically or convert manually before re-run")
    }
  }
  df <- psmelt(physeq)
  df <- df[, c("Abundance", variable, rank_names(physeq)[1:8])]
  df[, 2] <- gsub(" |\\(|\\)", "", df[, 2])
  df[, 2] <- as.factor(df[, 2])
  dir.create(output)
  for (lvl in levels(df[, 2])) {
    write.table(unique(df[which(df[, 2] == lvl), -2]), file = paste0(output, 
                                                                     "/", lvl, "taxonomy.txt"), sep = "\t", row.names = F, 
                col.names = F, na = "", quote = F)
  }
  krona_args <- paste(output, "/", levels(df[, 2]), "taxonomy.txt,", 
                      levels(df[, 2]), sep = "", collapse = " ")
  output <- paste(output, ".html", sep = "")
  system(paste("ktImportText", krona_args, "-o", output, sep = " "))
  browseURL(output)
}