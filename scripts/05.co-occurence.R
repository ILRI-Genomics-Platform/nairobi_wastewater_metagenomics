#Installing packages 
setwd("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/")
lib="/home/gkibet/R/x86_64-pc-linux-gnu-library/4.3"
source("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/scripts/functions.R")
requiredCRANPackages= c("BiocManager","dplyr", "ggplot2", "readr", "openxlsx","stringr", "xml2", "tidyverse",
                    "magrittr","ggVennDiagram","devtools","janitor","RColorBrewer","BBmisc","scales","Hmisc",
                    "igraph", "tidygraph", "graphlayouts", "ggraph","psych","RCy3")#,"RCytoscape"
requiredBIOCPackages=c("RCy3") #"RCy3",RCytoscape
#  CRAN packages will be installed if they are not yet installed.
installNloadCRANpackages(requiredPackages = requiredCRANPackages, lib = lib)
# BiocManager and GitHub packages will be installed using the functions below if they are not yet installed. Uncomment to install.
# installNloadBIOCpackages(requiredPackages = requiredBIOCPackages, lib = lib)
# installNloadGitHubpackages(requiredPackages = requiredGITHUBPackages, lib = lib)
# getwd() #Get working directory

shortDate <- gsub("-","",base::Sys.Date())
shortDate="20241002RAM"
#shortDate
#shortDate="20221116"


##Loading RGI data to R
taxa="family" #class|order|family|genera|species|
taxaAbundanceCutOff=2000
phage="phageGenera" #phageFamily|phageGenera|phageSpecies
phageAbundanceCutOff=50
rgiTerm="AMRGeneFamily" #AMRGeneFamily|AROTerm
coverageCutOff=50 #Alternative = log(AveragePercentCoverage)
coverage=paste(coverageCutOff,"PercCoverage",sep = "") #50PercCoverage|AnyPercCoverage|80PercCoverage
rgiAbundanceCutOff=20
rgiData_old <- read.csv("rgi/args/20240723_all_ARG-Metadata-10coverage.csv", header = T, sep = "\t")
rgiData000 <- read.csv("rgi/args/20240723_all_ARG-Metadata-50coverage.csv", header = T, sep = "\t")
rgiData00 <- read.csv("rgi/args/20240723_all_ARG-Metadata.csv", header = T, sep = "\t")
nameRename <- c(Name = "sampleID", sampleID = "Name")
# Renaming columns
rgiData <- rgiData00 %>% rename_with(~str_remove_all(., '\\.')) %>% filter(!is.na(.[[rgiTerm]])) %>%
  subset(AveragePercentCoverage >= coverageCutOff) #subset(log(AveragePercentCoverage) >= coverageCutOff)
hist(log(rgiData00$Average.Percent.Coverage))
hist(log(rgiData$AveragePercentCoverage))
# Check the data for those with 'NAs' in 'Name'
# rgiData001 <- rgiData[is.na(rgiData00$Name),]
# rgiData001 <- rgiData[is.na(rgiData$AMRGeneFamily),]

# Selecting RGI gene families
rgicleanDataRaw00 <- rgiData %>% select("Name",rgiTerm[1],"AllMappedReads","trimmedtotal_reads") %>% 
  group_by(Name,.[rgiTerm],trimmedtotal_reads) %>% 
  summarise_at(vars(AllMappedReads), list(AllMappedReads = sum)) %>% as.data.frame()
# Extract RGI data into a matrix
rgicleanDataRaw <- rgicleanDataRaw00 %>% select(-trimmedtotal_reads) %>%
  reshape(., idvar = rgiTerm[1], timevar="Name",v.names = "AllMappedReads",direction = "wide") %>%
  rename_with(~str_remove(., 'AllMappedReads.')) %>% rename(Name = rgiTerm[1]) %>% as.data.frame() %>%
  rowwise() %>% mutate(MeanAbundance=rowMeans(across(where(is.numeric)),na.rm = TRUE), .after = Name) %>%
  subset(MeanAbundance >= rgiAbundanceCutOff) %>% select(-MeanAbundance) #rgiAbundanceCutOff
cat("Original",rgiTerm[1],"count was",length(unique(rgicleanDataRaw00$AMRGeneFamily)),
    "\nAfter setting",rgiAbundanceCutOff,"minimum mean abundance cutoff we have",length(unique(rgicleanDataRaw$Name)),"\n")

# Extract Sample QC data into a matrix
rgiSampleQCdata <- rgicleanDataRaw00 %>% select(-rgiTerm[1]) %>% group_by(Name,trimmedtotal_reads) %>% 
  summarise_at(vars(AllMappedReads), list(AllMappedReads = sum)) %>% t() %>% as.data.frame() %>% 
  row_to_names(1) %>% mutate_if(is.character,as.numeric) %>% rownames_to_column(var = "Name")
rgiSampleQCdata00 <- rgiSampleQCdata %>% filter(Name == "trimmedtotal_reads") %>% select(-Name) %>%
  pivot_longer(cols = everything(), names_to = "sampleIDs", values_to = "trimmedtotal_reads")
# Generate a Relative abundance matrix of mapped RGI reads against trimmed_reads per sample
rgicleanDataRAb <- rgicleanDataRaw %>% 
  pivot_longer(cols = -c(1), names_to = "sampleIDs", values_to = "abundance") %>% 
  left_join(., rgiSampleQCdata00, by = c('sampleIDs')) %>% 
  mutate(rAbundance = as.numeric(abundance)/as.numeric(trimmedtotal_reads)) %>%
  select(Name, sampleIDs, rAbundance) %>% pivot_wider(names_from = sampleIDs, values_from = rAbundance)
# rgicleanDataRAb0 <- rgicleanDataRaw00 %>% 
#   mutate(across(-c(1,2,3),.fns=~./trimmedtotal_reads)) %>% select(-trimmedtotal_reads) %>% 
#   reshape(., idvar = rgiTerm[1], timevar="Name",v.names = "AllMappedReads",direction = "wide") %>%
#   rename_with(~str_remove(., 'AllMappedReads.')) %>% rename(Name = rgiTerm[1])
rgicleanDataRAMatrix <- rgicleanDataRAb %>% t() %>% as.data.frame() %>% row_to_names(1) %>%
  mutate_if(is.character,as.numeric) %>% rownames_to_column(var = "sampleID") %>% 
  mutate_at(c(-1), ~(scale(.,center=FALSE, scale=sd(., na.rm = TRUE))) %>% as.vector) %>% t() %>% as.data.frame() %>% 
  row_to_names(1) %>% mutate_if(is.character,as.numeric) %>% rownames_to_column(var = "Name")
# Generate a Relative abundance matrix of mapped RGI reads against sum of allmappedreads per sample
rgicleanDataPerc <- rgicleanDataRaw %>% as.data.frame() %>% mutate(across(-1, ~./sum(as.numeric(.),na.rm = TRUE)))

#RGI Metadata
rgiData %>% select(rgiTerm[1],"ResistanceMechanism") %>% distinct() %>% 
  group_by(.[rgiTerm]) %>% slice(1) -> rgiMetadataRaw
write.table(rgiMetadataRaw, file = paste("./rgi/plots/co-occurence/",shortDate,"_",rgiTerm,coverage,rgiAbundanceCutOff,"abundanceCutOff_","Unique.csv", sep = ""),
            row.names = F, col.names= TRUE, sep = "\t")
names(rgiMetadataRaw)[names(rgiMetadataRaw) == rgiTerm] <- "Name"
RGIData <- read.csv("plots/co-occurence/rgi/AMRGeneFamilyUnique_edited.csv", header = T, sep = "\t")
rgiMetadata <- left_join(rgiMetadataRaw,RGIData) %>% mutate(Name = ifelse(!is.na(shortName), shortName, Name)) %>%
  select("Name") %>% mutate(Group = "ARG") #

#Filtering rgicleanData
rgicleanPlotData=rgicleanDataRAMatrix #rgicleanDataPerc|rgicleanDataRAMatrix|rgicleanDataRAb
rgiINdata="rgicleanDataRAMatrix"
max_ARGs = 1000 #NULL|25|50|100|150|1000
minARG_freq = 50 #NULL|25|50|70|100|150
# Count the number of occurrences of ARO.Term in different samples.
rgicleanPlotData00 <- rgicleanPlotData %>% mutate(ARO.Term.Freq = rowSums(!is.na(.[,-1])), .after = Name)
# Calculate the mean of each ARO.Term across different sample.
rgicleanDataMean <- data.frame(Name=rgicleanPlotData00[,1], ARO.Term.Freq=rgicleanPlotData00[,2], 
                               Means=as.numeric(rowMeans(rgicleanPlotData00[,-1:-2], na.rm=T)))
# rgicleanDataMean00 <-rgicleanDataMean
hist(rgicleanDataMean$ARO.Term.Freq)
hist(log(rgicleanDataMean$Means))

# Filter ARO.Terms based on their frequency of occurrences: ARO.Term.Freq > minARG_freq 
rgicleanDataMeanFiltered <- rgicleanDataMean %>% subset(ARO.Term.Freq > minARG_freq) %>% 
  slice_max(order_by = as.numeric(Means), n = max_ARGs) 
rgicleanData00 <- semi_join(rgicleanPlotData,rgicleanDataMeanFiltered, by = "Name")
Scatter_Matrix <- pairs.panels(rgicleanData00[, c(2:21)], main = paste("Scatter Plot Matrix for rgi Dataset",sep=""))
pdf(paste("./rgi/plots/co-occurence/",shortDate,"_",rgiTerm,"_",rgiINdata,"_",coverage,rgiAbundanceCutOff,
               "abundanceCutOff","_SPLOM.pdf", sep = ""), width = 20, height = 20, fonts = "Helvetica", pointsize = 20)
pairs.panels(rgicleanData00[, c(2:21)], main = "Scatter Plot Matrix for rgi Dataset")
dev.off()

rgi_hist <- rgicleanData00 %>% pivot_longer(cols = -Name) %>% ggplot(aes(value)) +
  facet_wrap(~ Name, scales = "free") + 
  geom_histogram(bins = 30, fill = "lightblue",color = "black") + 
  geom_density(color = "black")
ggsave(paste("./rgi/plots/co-occurence/",shortDate,"_",rgiTerm,"_",rgiINdata,"_",coverage,rgiAbundanceCutOff,"abundanceCutOff","_scatterPlot.pdf", sep = ""),
       plot=rgi_hist, width = 20, height = 20, units = "in") 

# 
rgicleanData01 <- merge(rgicleanData00,select(RGIData,Name,shortName), by.x = "Name", by.y = "Name", all.x = TRUE) %>% 
  relocate(shortName, .after = Name) #%>% mutate(shortName = coalesce(shortName, Name))
# rgicleanData001 <- rgicleanData01 %>% group_by(Name) %>% replace(is.na(.$shortName), .$Name)
# rgicleanData01$shortName <- rowwise(replace(rgicleanData01$shortName, is.na(rgicleanData01$shortName), rgicleanData01$Name))
rgicleanData <- rgicleanData01 %>% relocate(shortName) %>% select(-c(shortName)) %>% rename_at(1,~"Name")
#hist(rgicleanDataMean$Means)
#rowMeans(rgicleanData[,-1], na.rm=T) %>% [order()]

# Extract ARO Terms related to selected AMRGeneFamilies:
selected_AMRGeneFamilies <- rgicleanData$Name
selected_AROTermsdf <- rgiData %>% filter(AMRGeneFamily %in% selected_AMRGeneFamilies) %>%
  select(AROTerm,AMRGeneFamily,DrugClass,ResistanceMechanism) %>% distinct()
write.table(selected_AROTermsdf, file = paste("./rgi/plots/co-occurence/",shortDate,"_",rgiTerm,coverage,rgiAbundanceCutOff,"abundanceCutOff_","selected_AROTerms.csv", sep = ""),
            row.names = F, col.names= TRUE, sep = "\t")
selected_AROTerms <- selected_AROTermsdf %>% select(AROTerm)

##Loading Taxa Data
taxaData0 <- read.csv(paste("rgi/abundance/bacteria/20241002_bacteriaMergedCladereads.csv", sep = ""), header = T, sep = "\t") %>% 
  select(-ends_with(".x")) %>% rename_with(~str_remove(., '.y')) %>% rename_with(~str_remove(., '.kraken2'))
taxaData00 <- read.csv(paste("rgi/abundance/bacteria/20240723_bacteriaMergedCladereads.csv", sep = ""), header = T, sep = "\t") %>% 
  select(-ends_with(".x")) %>% rename_with(~str_remove(., '.y')) %>% rename_with(~str_remove(., '.kraken2'))
drop_columns <- c("TaxRank","TaxID","OVERVIEW","STAT","TaxLineage","COVG004","COVG005")
#c("TaxRank","TaxID","OVERVIEW","STAT","TaxLineage","COVG00024_S8","COVG00090_S9","COVG00156_S10")
drop_rows <- c("root","unclassified")
taxaData <- taxaData0 %>% subset(TaxRank == 'F')
taxacleanDataRaw0 <- taxaData[!(taxaData$Name %in% drop_rows),!(names(taxaData) %in% drop_columns)] %>%
  rowwise() %>% mutate(MeanAbundance=rowMeans(across(where(is.numeric)),na.rm = TRUE)) 
taxacleanDataRaw <- taxacleanDataRaw0 %>% subset(MeanAbundance >= taxaAbundanceCutOff) %>% select(-MeanAbundance)
hist(taxacleanDataRaw0$MeanAbundance)
hist(log(taxacleanDataRaw0$MeanAbundance))
cat("Original bacterial",taxa,"count was",length(unique(taxaData$Name)),
    "\nAfter setting",taxaAbundanceCutOff,"minimum mean abundance cutoff we have",length(unique(taxacleanDataRaw$Name)),"\n")

taxaSampleQCdata <- taxaData0[(taxaData0$Name %in% drop_rows),!(names(taxaData0) %in% drop_columns)] %>%
  mutate(Name = replace(Name, Name == "root", "root_bacteria"))
sampleQCdata <- bind_rows(taxaSampleQCdata,rgiSampleQCdata)

# Generate a Relative abundance matrix of mapped taxa reads against trimmed_reads per sample
taxacleanDataRAb <- taxacleanDataRaw %>% 
  pivot_longer(cols = -c(1), names_to = "sampleIDs", values_to = "abundance") %>% 
  left_join(., sampleQCdata %>% filter(.,  Name == "trimmedtotal_reads") %>% select(-Name) %>%
              pivot_longer(cols = everything(), names_to = "sampleIDs", values_to = "abundance"), by = c('sampleIDs')) %>% 
  mutate(abundance = as.numeric(abundance.x)/as.numeric(abundance.y)) %>%
  select(Name, sampleIDs, abundance) %>% pivot_wider(names_from = sampleIDs, values_from = abundance)
taxacleanDataRAMatrix <- taxacleanDataRAb %>% t() %>% as.data.frame() %>% row_to_names(1) %>% 
  mutate_if(is.character,as.numeric) %>% rownames_to_column(var = "sampleID") %>% 
  mutate_at(vars(-sampleID), ~scale(.,center=FALSE, scale=sd(., na.rm = TRUE)) %>% as.vector) %>% t() %>%
  as.data.frame() %>% row_to_names(1) %>% mutate_if(is.character,as.numeric) %>% rownames_to_column(var = "Name")

# Generate a Relative abundance matrix of mapped taxa reads against sum of abundance per sample
taxacleanDataPerc <- taxacleanDataRaw %>% mutate(across(-1, ~./sum(as.numeric(.),na.rm = TRUE)))

#Taxa Metadata
taxaData %>% select("Name","TaxRank") %>% distinct() %>% mutate(Group = "Taxa") -> taxaMetadata

#Filtering taxacleanData
taxacleanPlotData=taxacleanDataRAMatrix #taxacleanDataPerc|taxacleanDataRAMatrix|taxacleanDataRAb
taxaINdata="taxacleanDataRAMatrix"
max_Taxa = 1000 #NULL|25|50|100|150|1000
minTaxa_freq = 50 #NULL|25|50|100|150
taxacleanPlotData00 <- taxacleanPlotData %>% mutate(taxa.Freq = rowSums(!is.na(.[,-1])), .after = Name)
taxacleanDataMean <- data.frame(Name=taxacleanPlotData00[,1], taxa.Freq=taxacleanPlotData00[,2], 
                               Means=as.numeric(rowMeans(taxacleanPlotData00[,-1:-2], na.rm=T)))
hist(taxacleanDataMean$taxa.Freq)
hist(log(taxacleanDataMean$Means))
taxacleanDataMeansFiltered <- taxacleanDataMean %>% subset(taxa.Freq > minTaxa_freq) %>% 
  slice_max(order_by = as.numeric(Means), n = max_Taxa) 
taxacleanData <- semi_join(taxacleanPlotData,taxacleanDataMeansFiltered,by = "Name")
Scatter_Matrix <- pairs.panels(taxacleanData[, c(2:10)], main = "Scatter Plot Matrix for rgi Dataset")
pdf(paste("./rgi/plots/co-occurence/",shortDate,"_","bacteria_",taxaINdata,"_SPLOM.pdf", sep = ""), 
    width = 20, height = 20, fonts = "Helvetica", pointsize = 20)
pairs.panels(taxacleanData[, c(2:21)], main = "Scatter Plot Matrix for rgi Dataset")
dev.off()

taxa_hist <- taxacleanData %>% pivot_longer(cols = -Name) %>% ggplot(aes(value)) +
  facet_wrap(~ Name, scales = "free") + 
  geom_histogram(bins = 30, fill = "lightblue",color = "black") + 
  geom_density(color = "black")
ggsave(paste("./rgi/plots/co-occurence/",shortDate,"_","bacteria_",taxaINdata,"_scatterPlot.pdf", sep = ""),
       plot=taxa_hist, width = 20, height = 20, units = "in") 

##Loading phage Data
phageData0 <- read.csv("rgi/abundance/phages/20241002_phageMergedCladereads.csv", header = T, sep = "\t") %>%
  select(-ends_with(".x")) %>% rename_with(~str_remove(., '.y')) %>% rename_with(~str_remove(., '.kraken2'))
# Filtering Genus data and samples
drop_columns <- c("TaxRank","TaxID","OVERVIEW","STAT","TaxLineage","COVG004","COVG005")
#c("TaxRank","TaxID","OVERVIEW","STAT","TaxLineage","COVG00024_S8","COVG00090_S9","COVG00156_S10")
drop_rows <- c("root","unclassified")
phageData <- phageData0 %>% subset(TaxRank == 'G')
#phageGeneraraw = phageTaxaraw
phageDataraw <- phageData[!(phageData$Name %in% drop_rows),!(names(phageData) %in% drop_columns)] %>%
  rowwise() %>% mutate(MeanAbundance=rowMeans(across(where(is.numeric)),na.rm = TRUE)) %>%
  subset(MeanAbundance >= phageAbundanceCutOff) %>% select(-MeanAbundance)
cat("Original phage",phage,"count was",length(unique(phageData$Name)),
    "\nAfter setting",phageAbundanceCutOff,"minimum mean abundance cutoff we have",length(unique(phageDataraw$Name)),"\n")

phageSampleQCdata <- phageData0[(phageData0$Name %in% drop_rows),!(names(phageData0) %in% drop_columns)] %>%
  mutate(Name = replace(Name, Name == "root", "root_phage"))
sampleMergedQCdata <- bind_rows(sampleQCdata,phageSampleQCdata) %>% distinct()

# Generate a Relative abundance matrix of mapped taxa reads against trimmed_reads per sample
phagecleanDataRAb <- phageDataraw %>% 
  pivot_longer(cols = -c(1), names_to = "sampleIDs", values_to = "abundance") %>% 
  left_join(., sampleQCdata %>% filter(.,  Name == "trimmedtotal_reads") %>% select(-Name) %>%
              pivot_longer(cols = everything(), names_to = "sampleIDs", values_to = "abundance"),
            by = c('sampleIDs')) %>% 
  mutate(abundance = as.numeric(abundance.x)/as.numeric(abundance.y)) %>%
  select(Name, sampleIDs, abundance) %>% pivot_wider(names_from = sampleIDs, values_from = abundance)
phagecleanDataRAMatrix <- phagecleanDataRAb %>% t() %>% as.data.frame() %>% row_to_names(1) %>%
  mutate_if(is.character,as.numeric) %>% rownames_to_column(var = "sampleID") %>% 
  mutate_at(c(-1), ~(scale(.,center=FALSE, scale=sd(., na.rm = TRUE))) %>% as.vector) %>% t() %>% as.data.frame() %>% 
  row_to_names(1) %>% mutate_if(is.character,as.numeric) %>% rownames_to_column(var = "Name")

# Generate a Relative abundance matrix of mapped taxa reads against sum of abundance per sample
phagecleanTaxaPerc <- phageDataraw %>% mutate(across(-1, ~./sum(as.numeric(.),na.rm = TRUE)))

#Phage Metadata
phageData %>% select("Name","TaxRank") %>% distinct() %>% mutate(Group = "phage") -> phageMetadata

#Filtering phageGeneraPerc
phagecleanPlotData=phagecleanDataRAMatrix #phagecleanTaxaPerc|phagecleanDataRAMatrix|phagecleanDataRAb
phageINdata="phagecleanDataRAMatrix"
max_Phage = 1000 #NULL|25|50|100|150|1000
minPhage_freq = 50 #NULL|25|50|100|150
# max_Phage = 10 #
phagecleanPlotData00 <- phagecleanPlotData %>% mutate(phage.Freq = rowSums(!is.na(.[,-1])), .after = Name)
phagecleanDataMean <- data.frame(Name=phagecleanPlotData00[,1], phage.Freq=phagecleanPlotData00[,2], 
                                Means=as.numeric(rowMeans(phagecleanPlotData00[,-1:-2], na.rm=T)))
phagecleanDataMeansFiltered <- phagecleanDataMean %>% subset(phage.Freq > minPhage_freq) %>% 
  slice_max(order_by = as.numeric(Means), n = max_Phage) 
phageCleanData <- semi_join(phagecleanPlotData,phagecleanDataMeansFiltered,by = "Name")
Scatter_Matrix <- pairs.panels(phageCleanData[, c(2:10)], main = "Scatter Plot Matrix for rgi Dataset")
pdf(paste("./rgi/plots/co-occurence/",shortDate,"_","bacteria_",phageINdata,"_SPLOM.pdf", sep = ""), 
    width = 20, height = 20, fonts = "Helvetica", pointsize = 20)
pairs.panels(taxacleanData[, c(2:21)], main = "Scatter Plot Matrix for rgi Dataset")
dev.off()

phage_hist <- phageCleanData %>% pivot_longer(cols = -Name) %>% ggplot(aes(value)) +
  facet_wrap(~ Name, scales = "free") +
  geom_histogram(bins = 30, fill = "lightblue",color = "black") +
  geom_density(color = "black")
ggsave(paste("./rgi/plots/co-occurence/",shortDate,"_",phageINdata,"_scatterPlot.pdf", sep = ""),
       plot=phage_hist, width = 20, height = 20, units = "in",limitsize = FALSE)

# Merging Taxa and RGI Data
taxa.Data=taxacleanData
rgi.Data=rgicleanData
phage.Data = phageCleanData
setdiff(colnames(taxa.Data),colnames(rgi.Data))
setdiff(colnames(taxa.Data),colnames(phage.Data))
cleanData00 <- t(bind_rows(taxa.Data,rgi.Data)) %>% row_to_names(1) 
cleanData01 <- t(bind_rows(taxa.Data,rgi.Data,phage.Data)) %>% row_to_names(1)

cleanData00[is.na(cleanData00)] <- 0
cleanData01[is.na(cleanData01)] <- 0

#Merging Taxa and RGI metadata and filtering
metadataRaw <- bind_rows(taxaMetadata, rgiMetadata, phageMetadata)
#bind_rows(taxa.Data,rgi.Data) -> cleanDataJ
cleanMetadata <- semi_join(metadataRaw, bind_rows(taxa.Data,rgi.Data, phage.Data), by = "Name") %>% 
  as.data.frame() %>% mutate_at(2, ~replace_na(.,rgiTerm))
#setdiff(cleanMetadata$Name,cleanDataJ$Name)

#Plot Data
# cleanData=cleanData00 %>% as.data.frame() %>% mutate_if(is.character,as.numeric)
cleanData=cleanData00 %>% as.data.frame() %>% mutate_if(is.character,as.numeric)
if (rgiINdata == "rgicleanDataRAMatrix") {
  INdata="_cleanDataRAMatrixTRP_"
} else if (rgiINdata == "rgicleanDataRAb") {
  INdata="_cleanDataRAbTRP_"
} else {
  INdata="_cleanDataTaxaPercTRP_"
}
pairs.panels(cleanData[, c(2:21)], main = "Scatter Plot Matrix for rgi-bacteria-phage Dataset")
pdf(paste("./rgi/plots/co-occurence/",shortDate,INdata,"bacteria_",taxaINdata,rgiINdata,phageINdata,"_SPLOM.pdf", sep = ""), 
    width = 100, height = 100, fonts = "Helvetica", pointsize = 20)
minVal=length(colnames(cleanData))-100
maxVal=length(colnames(cleanData))
pairs.panels(cleanData[, c(minVal:maxVal)], main = "Scatter Plot Matrix for rgi Dataset")
dev.off()

#Spearman's Co-rrelation Coefficients - Using dplyr::cor() function
rho_coeff <- cor(cleanData, method = "spearman")
rho_coeffClean <- ifelse(rho_coeff<=0.6,0,rho_coeff)
rho_coeffdf <- setNames(as.data.frame(as.table(rho_coeff)), c('Var1','Var2','weight'))
rhoNet <- graph_from_adjacency_matrix(rho_coeffClean, weighted=T, mode="undirected", diag=F)
#plot(rhoNet)

#Spearman's Co-rrelation Coefficients - Using hmisc::rcorr() function
rho_coeff.rcorr <- rcorr(as.matrix(cleanData), type = "spearman")
rho_coeff.pvalue <- round(rho_coeff.rcorr[["P"]],6)
rho_coeff.pvalue.df <- setNames(as.data.frame(as.table(rho_coeff.pvalue)), c('Var1','Var2','P'))
rho_coeff.rvalue <- rho_coeff.rcorr[["r"]]
rho_coeff.rvalue.df <- setNames(as.data.frame(as.table(rho_coeff.rvalue)), c('Var1','Var2','weight'))
rho_coeff.rcorr.df <- full_join(rho_coeff.rvalue.df, rho_coeff.pvalue.df, by = c("Var1","Var2")) #%>% setNames(c("Var1","Var2","weight","P")) 
rho_coeff.rcorr.df.filtered <- filter(rho_coeff.rcorr.df, weight > 0.75 & P < 0.05) 
rho_coeff.rcorr.df.Clean <- rho_coeff.rcorr.df.filtered %>% select(-P) #%>% 
#   merge(.,select(cleanMetadata,Name,Group), by.x = "Var1", by.y = "Name", all.x = TRUE) %>%
#   merge(.,select(cleanMetadata,Name,Group), by.x = "Var2", by.y = "Name", all.x = TRUE) %>%
#   subset(Group.x != Group.y) %>% select(-c("Group.x","Group.y"))
rho_coeff.rcorr.clean <- rho_coeff.rcorr.df.Clean %>% 
  reshape(., idvar = "Var1", timevar = "Var2", direction = "wide") %>%
  set_rownames(.$Var1) %>% rename_with(~str_remove_all(., 'weight.')) %>% 
  select(2:length(colnames(.))) %>% mutate_all(.,~ifelse(is.na(.), 0,.)) %>% as.matrix()
rhoPvalNet <- graph_from_adjacency_matrix(rho_coeff.rcorr.clean, weighted=T, mode="undirected", diag=F)
#plot(rhoPvalNet)

#Co-occurence Network graph
network = rhoPvalNet #rhoNet|rhoPvalNet
argnetclean = "RhoPval" #Rho|Pval|RhoPval
netClean = paste("_",coverage,argnetclean,"Filtered",sep="") #Rho|Pval|RhoPval
#plot(network)

#Plot a rudimentary network - jpeg
jpeg(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
           rgiTerm,"Top",max_ARGs,netClean,"_clusterbetweeness0",".jpeg", sep = ""),
     width = 480, height = 480,units = "mm", res=200)
dev.set(which = 2)
plot(network)
dev.copy(which = 4)
dev.off()
#Plot a rudimentary network - pdf
pdf(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
          rgiTerm,"Top",max_ARGs,netClean,"_clusterbetweenesso",".pdf", sep = ""),
    width = 480, height = 480, fonts = "Helvetica", pointsize = 400)
plot(network)
dev.off()


#Write network data to file
write.table(rho_coeff, file = paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",
                                    phage,"Top",max_Phage,"_",rgiTerm,"Top",max_ARGs,netClean,"-rhoAll.csv", sep = ""),
            row.names = TRUE, col.names= TRUE, sep = "\t")
write.table(rho_coeffClean, file = paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",
                                         phage,"Top",max_Phage,"_",rgiTerm,"Top",max_ARGs,netClean,"-rhoClean.csv", sep = ""),
            row.names = TRUE, col.names= TRUE, sep = "\t")
write.table(rho_coeffdf, file = paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",
                                      phage,"Top",max_Phage,"_",rgiTerm,"Top",max_ARGs,netClean,"-rhoAllTable.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.table(rho_coeff.rcorr.clean, file = paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",
                                                phage,"Top",max_Phage,"_",rgiTerm,"Top",max_ARGs,netClean,"-rhoPvalClean.csv", sep = ""),
            row.names = TRUE, col.names= TRUE, sep = "\t")
write.table(rho_coeff.pvalue, file = paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",
                                           phage,"Top",max_Phage,"_",rgiTerm,"Top",max_ARGs,netClean,"-PvalAllTable.csv", sep = ""),
            row.names = TRUE, col.names= TRUE, sep = "\t")

#ggplot option
#ggraph(network, layout = 'auto') +
#  geom_edge_link(alpha = 0.25) +
#  geom_node_point(color="lightblue")

# Check cluster Edge betweenness - classic Girvan & Newman betweenness clustering method
# Community structure detection algorithms try to find dense subgraphs within larger network graphs
cb <- cluster_edge_betweenness(network)
jpeg(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
           rgiTerm,"Top",max_ARGs,netClean,"_clusterbetweeness",".jpeg", sep = ""),
     width = 480, height = 480,units = "mm", res=200)
dev.set(which = 2)
plot(cb, y=network, vertex.label=NULL,  vertex.size=3)
dev.copy(which = 4)
dev.off()
pdf(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
          rgiTerm,"Top",max_ARGs,netClean,"_clusterbetweeness",".pdf", sep = ""),
     width = 480, height = 480, pointsize = 400)
plot(cb, y=network, vertex.label=NULL,  vertex.size=3)
dev.off()

# extract a cluster/community membership vector for further inspection with the membership() function:
head( membership(cb) )

#Node degree - Number of adjacent nodes to a node/vertex
dg <- degree(network)
hist(dg, breaks=30, col="lightblue", main ="Node Degree Distribution")

jpeg(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
           rgiTerm,"Top",max_ARGs,netClean,"_NodeDegree",".jpeg", sep = ""),
     width = 200, height = 200,units = "mm", res=200)
dev.set(which = 2)
plot(degree_distribution(network),type="h")
dev.copy(which = 4)
dev.off()
pdf(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
          rgiTerm,"Top",max_ARGs,netClean,"_NodeDegree",".pdf", sep = ""),
     width = 200, height = 200, pointsize = 150)
plot(degree_distribution(network),type="h")
dev.off()


#Centrality Analysis - Estimate how important a node/edge is for connectivity of the network
#Using Google PageRank
rho_rank <- page_rank(network)
v.size <- BBmisc::normalize(rho_rank$vector, range=c(2,20), method="range")
# View(as.data.frame(v.size))

jpeg(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",rgiTerm,"Top",
           max_ARGs,netClean,"_PageRankCA",".jpeg", sep = ""),
     width = 480, height = 480,units = "mm", res=200)
dev.set(which = 2)
plot(network, vertex.size=v.size, vertex.label=NA)
dev.copy(which = 4)
dev.off()
pdf(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
          rgiTerm,"Top",max_ARGs,netClean,"_PageRankCA",".pdf", sep = ""),
     width = 480, height = 480, pointsize = 300)
plot(network, vertex.size=v.size, vertex.label=NA)
dev.off()

#CA Using Node Degree (dg)
rho_NodeDegree <- degree(network)
v.size <- BBmisc::normalize(rho_NodeDegree, range=c(2,20), method="range")
jpeg(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
           rgiTerm,"Top",max_ARGs,netClean,"_NodeDegreeCA",".jpeg", sep = ""),
     width = 480, height = 480,units = "mm", res=200)
dev.set(which = 2)
plot(network, vertex.size=v.size, vertex.label=NA)
dev.copy(which = 4)
dev.off()
pdf(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
          rgiTerm,"Top",max_ARGs,netClean,"_NodeDegreeCA",".pdf", sep = ""),
     width = 480, height = 480, pointsize = 300)
plot(network, vertex.size=v.size, vertex.label=NA)
dev.off()

#Centrality score based on Edge betweeness
rho_clusterBetweeness <- betweenness(network)
v.size <- BBmisc::normalize(rho_clusterBetweeness, range=c(2,20), method="range")

jpeg(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
           rgiTerm,"Top",max_ARGs,netClean,"_ClusterBetweenessCA",".jpeg", sep = ""),
     width = 480, height = 480,units = "mm", res=200)
dev.set(which = 2)
plot(network, vertex.size=v.size)
dev.copy(which = 4)
dev.off()
pdf(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
          rgiTerm,"Top",max_ARGs,netClean,"_ClusterBetweenessCA",".pdf", sep = ""),
     width = 480, height = 480, pointsize = 400)
plot(network, vertex.size=v.size)
dev.off()

genenet.nodes <- data.frame(subset(cleanMetadata,Name %in% V(network)$name))
genenet.edges <- data.frame(igraph::as_edgelist(network))
genenet.edges$Weight <- igraph::edge_attr(network)[[1]]
names(genenet.edges) <- c("source","target","interaction")
genenet.edges$source <- as.character(genenet.edges$source)
genenet.edges$target <- as.character(genenet.edges$target)
names(genenet.nodes) <- c("id","TaxRank","Group")
genet.raw <- graph_from_data_frame(d=genenet.edges,vertices = genenet.nodes, directed = F)
# Can also use RCy3 -- see the below
# Identification of all nodes with less than 2 edges
verticesToRemove <- V(genet.raw)[degree(genet.raw) < 2]
# These edges are removed from the graph
genet <- delete.vertices(genet.raw, verticesToRemove) 
#Vertex size
v.size <- BBmisc::normalize(rho_clusterBetweeness, range=c(2,20), method="range")
rho_clusterBetweeness <- betweenness(genet)
# Assign colors to nodes (search term blue, others orange)
V(genet)$color <- ifelse(V(genet)$Group == "Taxa", 'cornflowerblue',
                         ifelse(V(genet)$Group == "ARG", 'darkorange',
                                ifelse(V(genet)$Group == "phage",'darkolivegreen3', 
                                       'darkgray')))
# Set node lebel colours:
V(genet)$color_l <- ifelse(V(genet)$Group == "Taxa", 'blue4',
                           ifelse(V(genet)$Group == "ARG", 'brown4',
                                  ifelse(V(genet)$Group == "phage", 'darkslategrey',
                                         'darkgray')))
# Set edge colors
E(genet)$color <- adjustcolor("darkorchid4", alpha.f = .5)
# scale significance between 1 and 10 for edge width
ecount(genet)
#E(genet)$width <- scales::rescale(E(genet)$interaction, to = c(1, 10))

# Set edges with radius
E(genet)$curved <- 0.30 
# Size the nodes by their degree of networking (scaled between 5 and 15)
#V(genet)$size <- scales::rescale(log(degree(genet)), to = c(5, 15))

# Define the frame and spacing for the plot
par(mai=c(0,0,1,0)) 

#ceb <- cluster_edge_betweenness(genet) 
#dendPlot(ceb,mode="hclust")

jpeg(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
           rgiTerm,"Top",max_ARGs,netClean,"_CBCA_edited","layout_with_lgl01.jpeg", sep = ""),
     width = 480, height = 480,units = "mm", res=300)
dev.set(which = 2)
plot(genet, 
     layout = layout_with_lgl, # Force Directed Layout 
     main = paste(rgiTerm,"-",taxa, 'Co-occurence Network'),
     vertex.label.family = "sans",
     vertex.label.cex = 2,
     vertex.shape = "circle",
     vertex.label.dist = 0.5,          # Labels of the nodes moved slightly
     vertex.frame.color = adjustcolor("Black", alpha.f = .5),
     vertex.label.color = V(genet)$color_l,     # Color of node names
     vertex.label = V(genet)$name,      # node names #V(genet)$name|c("")
     vertex.size = v.size,
     vertex.color = adjustcolor(V(genet)$color, alpha.f = .7),)
dev.copy(which = 4)
dev.off()

pdf(paste("./rgi/plots/co-occurence/",shortDate,INdata,taxa,"TaxaTop",max_Taxa,"_",phage,"Top",max_Phage,"_",
          rgiTerm,"Top",max_ARGs,netClean,"_CBCA_edited","layout_with_lgl01.pdf", sep = ""),
     width = 480, height = 480, pointsize = 300)
plot(genet, 
     layout = layout_with_lgl, # Force Directed Layout 
     main = paste(rgiTerm,"-",taxa, 'Co-occurence Network'),
     vertex.label.family = "sans",
     vertex.label.cex = 2,
     vertex.shape = "circle",
     vertex.label.dist = 0.5,          # Labels of the nodes moved slightly
     vertex.frame.color = adjustcolor("Black", alpha.f = .5),
     vertex.label.color = V(genet)$color_l,     # Color of node names
     vertex.label = V(genet)$name,  # node names #V(genet)$name|c("")
     vertex.size = v.size,
     vertex.color = adjustcolor(V(genet)$color, alpha.f = .7),)
dev.off()

#  RCy3
# Send network to Cytoscape using RCy3
# Open a new connection and delete any existing windows/networks in Cy
cy <- CytoscapeConnection()
deleteAllWindows(cy)

ug <- cyPlot(genenet.nodes,genenet.edges)

