#Installation and setup
setwd("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/")
# setwd("/home/gkibet/bioinformatics/github/metagenomics/data/20230120_UrbanZooProj_EF_NextSeqHT/visualization/")
lib="/home/gkibet/R/x86_64-pc-linux-gnu-library/4.3"
source("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/scripts/functions.R")

## Parsing command line argument
args = commandArgs(trailingOnly=TRUE)
ncores0 = as.numeric(args[1])
print(ncores0)
print(typeof(ncores0))
ncores=ncores0 - 2 #1|16|15
# ncores=18
cat("\n\tThe number of cores available to R:",ncores)

# The list of all packages to be loaded.
requiredCRANPackages= c("BiocManager", 'broom', 'bslib', 'cli', 'FactoMineR', 'formatR', 'highr', 'foreach',
                        'htmlwidgets', 'httpuv', 'isoband', 'tidyverse', 'pkgdown', 'purrr',"openxlsx", 'packrat',
                        'rmarkdown', 'rsconnect', 'shiny', 'shinyWidgets', 'stringi', 'svglite', 'xfun', 'yulab.utils', 
                        'pavian','maptools','doParallel','future','future.apply','parallel','doFuture') #,"remotes"
# The two below are used to install the non-CRAN packages from BiocManager and GitHub if not yet installed.
requiredBIOCPackages=c("Rsamtools")
requiredGITHUBPackages = c("fbreitwieser/pavian")

#  CRAN packages will be installed if they are not yet installed.
installNloadCRANpackages(requiredPackages = requiredCRANPackages, lib = lib)
# BiocManager and GitHub packages will be installed using the functions below if they are not yet installed. Uncomment to install.
# installNloadBIOCpackages(requiredPackages = requiredBIOCPackages, lib = lib)
# installNloadGitHubpackages(requiredPackages = requiredGITHUBPackages, lib = lib)
# install.packages("maptools", repos="http://R-Forge.R-project.org")

# To launch shiny app - a shiny App alternative to the code below.
# #Launching Pavian ShinyApp from R, type
# pavian::runApp(port=5000)

# Read metadata file
dataDate="20241002"
runs=c("run01","run02","run03","run04","run05","run08","run09")#"run00",
metadata <- read.xlsx("metadata/20230824_wastewater_samples_metadata.xlsx", colNames = T,
                      sheet = "sequence_data") %>% filter(.,Seqrun %in% runs)
# metadata <- read.csv("./metadata/20240721_sequence_metadata.csv", header = T, sep = ",")
#View(metadata)
shortDate <- gsub("-","",base::Sys.Date())
shortDate="20241002"
#shortDate

# Create input df for filter_Reports() - Must have two columns 'Name' and 'ReportFilePath'
# 'Name' column contains preferred names for the samples from which the report was generated
# 'ReportFilePath' column contains filepaths to the kreports matching the 'Name'
reportFiles0 <- metadata %>% select(Name,sampleID,ReportFilePath_kraken2,tool,Seqrun) %>% 
  mutate(kreportFilePath_Kraken2 = gsub('kreports/wastewater/kraken2', 'rgi/stats',
                               gsub('_S\\d*.kreport.txt','.ARG-Kraken2.kreport.txt', ReportFilePath_kraken2, 
                                    ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)))
reportFilesList <- list.files(path=paste(getwd(),"/rgi/stats",sep = ""),pattern = ".ARG-Kraken2.kreport.txt", full.names = T)
reportFiles <- reportFiles0 %>% filter(kreportFilePath_Kraken2 %in% reportFilesList)
reportFilesnon <- reportFiles0 %>% filter(!(kreportFilePath_Kraken2 %in% reportFilesList))
#reportFilesdf = reportFiles

#taxRanks to keep
filteredTaxRanks = c("D", "K", "P", "C", "O", "F", "G", "S")

#Retain Bacterial taxa only
outputDir="./rgi/abundance/"
filteredDomains = c("d_Archaea","d_Eukaryota","d_Viruses")
kreportList <- filter_Reports(reportFilesdf=reportFiles,filteredTaxRanks=filteredTaxRanks,filteredTaxa=filteredDomains,
                              tool="tool",sampleID="Name",ReportFilePath="kreportFilePath_Kraken2",nCores=ncores)
saveRDS(kreportList, file = paste(outputDir,"bacteria/",shortDate,"_bacteriaMergedReads.RData", sep = ""))
# kreportList <- readRDS(paste(outputDir,"bacteria/",shortDate,"_bacteriaMergedReads.RData", sep = ""))

#Merging samples
cat("\n\tDone filtering Bacterial taxonomic groups...\n",
    "\tWill proceed and write out the output to is: ",outputDir,"bacteria/",shortDate,"*")
bacteriaMergedTaxareads = data.frame()
bacteriaMergedCladereads = data.frame()
bacteriaMergedTaxareads <- merge_reports(kreportList, numeric_col = c("taxonReads"))
bacteriaMergedCladereads <- merge_reports(kreportList, numeric_col = c("cladeReads"))
#Writing output
write.table(bacteriaMergedTaxareads, file = paste(outputDir, "bacteria/",shortDate,"_bacteriaMergedTaxareads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(bacteriaMergedTaxareads, file = paste(outputDir,"bacteria/",shortDate,"_bacteriaMergedTaxareads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")
write.table(bacteriaMergedCladereads, file = paste(outputDir,"bacteria/",shortDate,"_bacteriaMergedCladereads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(bacteriaMergedCladereads, file = paste(outputDir,"/bacteria/",shortDate,"_bacteriaMergedCladereads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")


# #Retain Viral taxa only
Domains = c("d_Archaea","d_Eukaryota","d_Bacteria","d_Virus")

#Taxa to filter
filteredDomains = c("d_Archaea","d_Eukaryota","d_Bacteria") #d_Virus

#taxRanks to keep
filteredTaxRanks = c("D", "K", "P", "C", "O", "F", "G", "S")

#Filtering Viral reads
# kreportList <- filter_Reports(reportFilesdf=reportFiles,filteredTaxRanks=filteredTaxRanks,filteredTaxa=filteredDomains,
#                               tool="tool",sampleID="Name",ReportFilePath="kreportFilePath_Kraken2",nCores=ncores)
kreportList <- filter_Reports(reportFilesdf=reportFiles,filteredTaxRanks=filteredTaxRanks,filteredTaxa=filteredDomains,
                              tool="tool",sampleID="Name",ReportFilePath="kreportFilePath_Kraken2",nCores=ncores)
saveRDS(kreportList, file = paste(outputDir,"/viruses/",shortDate,"_viralMergedReads.RData", sep = ""))
# kreportList <- readRDS(paste(outputDir,"/viruses/",shortDate,"_viralMergedReads.RData", sep = ""))

#Merging samples
cat("\n\tDone filtering Viral taxonomic groups...",
    "\n\tWill proceed and write out the output to is: ",outputDir,"viruses/",shortDate,"*")
viralMergedTaxareads = data.frame()
viralMergedCladereads = data.frame()
viralMergedTaxareads <- merge_reports(kreportList, numeric_col = c("taxonReads"))
viralMergedCladereads <- merge_reports(kreportList, numeric_col = c("cladeReads"))
#Writing output
write.table(viralMergedTaxareads, file = paste(outputDir,"/viruses/",shortDate,"_viralMergedTaxareads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(viralMergedTaxareads, file = paste(outputDir,"/viruses/",shortDate,"_viralMergedTaxareads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")
write.table(viralMergedCladereads, file = paste(outputDir,"/viruses/",shortDate,"_viralMergedCladereads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(viralMergedCladereads, file = paste(outputDir,"/viruses/",shortDate,"_viralMergedCladereads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")

# Filtering Phage data minus Viral data 
viralTaxadf <- read.csv("/home/gkibet/bioinformatics/github/metagenomics/docs/ictv/viralTaxa_ictv.tsv", header = T, sep = "\t")
viralTaxa <- as.vector(viralTaxadf$pavian_name)
viralTaxaNoPhage <- filter(viralTaxadf, Group == "virus") %>% pull(pavian_name)
phageTaxa <- filter(viralTaxadf, Group == "phage") %>% pull(pavian_name)

Domains = c("d_Archaea","d_Eukaryota","d_Bacteria")
filteredDomains = c(Domains,viralTaxaNoPhage)

#taxRanks to keep
filteredTaxRanks = c("C","O", "F", "G", "S")

# #Filtering for phage reads
kreportList <- filter_Reports(reportFilesdf=reportFiles,filteredTaxRanks=filteredTaxRanks,filteredTaxa=filteredDomains,
                              tool="tool",sampleID="Name",ReportFilePath="kreportFilePath_Kraken2",nCores=ncores)
saveRDS(kreportList, file = paste(outputDir,"/phages/",shortDate,"_phageMergedReads.RData", sep = ""))
# kreportList <- readRDS("./plotdata/abundance/phages/20240502_phageMergedTaxareads.RData")

#Merging samples
cat("\n\tDone filtering Phages taxonomic groups...",
    "\n\tWill proceed and write out the output to is: ",outputDir,"phages/",shortDate,"*")
phageMergedTaxareads = data.frame()
phageMergedCladereads = data.frame()
phageMergedTaxareads <- merge_reports(kreportList, numeric_col = c("taxonReads"))
phageMergedCladereads <- merge_reports(kreportList, numeric_col = c("cladeReads"))
#Writing output
write.table(phageMergedTaxareads, file = paste(outputDir,"/phages/",shortDate,"_phageMergedTaxareads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(phageMergedTaxareads, file = paste(outputDir,"/phages/",shortDate,"_phageMergedTaxareads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")
write.table(phageMergedCladereads, file = paste(outputDir,"/phages/",shortDate,"_phageMergedCladereads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(phageMergedCladereads, file = paste(outputDir,"/phages/",shortDate,"_phageMergedCladereads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")


# Filtering Viral data minus phage data
# Taxa to filter
Domains = c("d_Archaea","d_Eukaryota","d_Bacteria") #d_Virus
# phageClass= 'c_Caudoviricetes'
filteredDomains = c(Domains,phageTaxa)

#taxRanks to keep
filteredTaxRanks = c("D", "K", "P", "C", "O", "F", "G", "S")

# #Filtering Viral reads
kreportList <- filter_Reports(reportFilesdf=reportFiles,filteredTaxRanks=filteredTaxRanks,filteredTaxa=filteredDomains,
                              tool="tool",sampleID="Name",ReportFilePath="kreportFilePath_Kraken2",nCores=ncores)
saveRDS(kreportList, file = paste(outputDir,"/viruses/",shortDate,"_viralNoPhageMergedReads.RData", sep = ""))
# kreportList <- readRDS("./plotdata/abundance/viruses/20240502_viralNoPhageMergedTaxareads.RData")

#Merging samples
cat("\n\tDone filtering Phages taxonomic groups...",
    "\n\tWill proceed and write out the output to is: ",outputDir,"viruses/",shortDate,"*")
viralNoPhageMergedTaxareads = data.frame()
viralNoPhageMergedCladereads = data.frame()
viralNoPhageMergedTaxareads <- merge_reports(kreportList, numeric_col = c("taxonReads"))
viralNoPhageMergedCladereads <- merge_reports(kreportList, numeric_col = c("cladeReads"))

#Writing output
write.table(viralNoPhageMergedTaxareads, file = paste(outputDir,"/viruses/",shortDate,"_viralNoPhageMergedTaxareads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(viralNoPhageMergedTaxareads, file = paste(outputDir,"/viruses/",shortDate,"_viralNoPhageMergedTaxareads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")
write.table(viralNoPhageMergedCladereads, file = paste(outputDir,"/viruses/",shortDate,"_viralNoPhageMergedCladereads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(viralNoPhageMergedCladereads, file = paste(outputDir,"/viruses/",shortDate,"_viralNoPhageMergedCladereads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")

# taxRanks to keep
filteredTaxRanks = c("D", "K", "P", "C", "O", "F", "G", "S")

# Retain allDomains
filteredDomains = c("")
kreportList <- filter_Reports(reportFilesdf=reportFiles,filteredTaxRanks=filteredTaxRanks,filteredTaxa=filteredDomains,
                              tool="tool",sampleID="sampleID",ReportFilePath="ReportFilePath",nCores=ncores)
# Merging samples
allDomainsMergedTaxareads = data.frame()
allDomainsMergedCladereads = data.frame()
allDomainsMergedTaxareads <- merge_reports(kreportList, numeric_col = c("taxonReads"))
allDomainsMergedCladereads <- merge_reports(kreportList, numeric_col = c("cladeReads"))
# Writing output
write.table(allDomainsMergedTaxareads, file = paste("./plotdata/abundance/allDomains/",shortDate,"_allDomainsMergedTaxareads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(allDomainsMergedTaxareads, file = paste("./plotdata/abundance/allDomains/",shortDate,"_allDomainsMergedTaxareads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")
write.table(allDomainsMergedCladereads, file = paste("./plotdata/abundance/allDomains/",shortDate,"_allDomainsMergedCladereads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(allDomainsMergedCladereads, file = paste("./plotdata/abundance/allDomains/",shortDate,"_allDomainsMergedCladereads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")

# #Alternative: Filtering target pathogens from ViralMergedData and BacteriallMergedData
# Reading the target pathogen TaxIDs
pathogenData <- read.xlsx("../../assets/pathogens/Pathogens_of_Public_Interest_WHO.xlsx", colNames = T, 
                          sheet = "pathogen_taxids") %>% lapply(., function(y) gsub('\\&#10;','',y)) %>%
  as.data.frame()
pathogenTaxIDs <- pathogenData %>% separate_rows(everything(), sep = ";") #%>% select(-Pathogens) 

# Read taxonomic classification data
viralMergedTaxareads <- read.csv("./plotdata/abundance/viruses/20230718_viralMergedTaxareads.csv", header = T, sep = "\t")
viralMergedCladereads <- read.csv("./plotdata/abundance/viruses/20230718_viralMergedCladereads.csv", header = T, sep = "\t")
bacteriaMergedTaxareads <- read.csv("./plotdata/abundance/bacteria/20230718_bacteriaMergedTaxareads.csv", header = T, sep = "\t")
bacteriaMergedCladereads <- read.csv("./plotdata/abundance/bacteria/20230718_bacteriaMergedCladereads.csv", header = T, sep = "\t")

# Merge Viral & Bacterial data:
MergedTaxareads <- bind_rows(bacteriaMergedTaxareads,viralMergedTaxareads) %>% distinct() %>% 
  mutate_at(c("TaxID"), as.character)
MergedCladereads <- bind_rows(bacteriaMergedCladereads,viralMergedCladereads) %>% distinct() %>%
  mutate_at(c("TaxID","OVERVIEW","STAT"), as.character) %>% group_by(pick(where(is.character))) %>% 
  summarise(across(matches("COVG"), sum), .groups = "drop") %>% select(colnames(bacteriaMergedCladereads))
MergedCladereads %>% group_by(TaxID) %>% filter(., TaxID %in% c(subset(.,duplicated(TaxID))$TaxID)) %>%
  arrange(Name) -> duplicateTaxa

# # Filtering out pathogens from merged data:
dropCols=c("OVERVIEW","STAT")
pathogenMergedTaxareads <- filter(MergedTaxareads, MergedTaxareads$TaxID %in% pathogenTaxIDs$TaxID) %>% 
  select(-dropCols) %>% left_join(., pathogenTaxIDs, by = join_by(TaxID)) %>% 
  relocate(Domain, .after = TaxID)
pathogenMergedCladereads <- filter(MergedCladereads, MergedCladereads$TaxID %in% pathogenTaxIDs$TaxID) %>% 
  select(-dropCols) %>% left_join(., pathogenTaxIDs, by = join_by(TaxID)) %>% 
  relocate(Domain, .after = TaxID)
#setdiff(pathogenTaxIDs$TaxID,pathogenMergedCladereads$TaxID)

write.table(pathogenMergedCladereads, file = paste(outputDir,"/pathogens/",shortDate,"_pathogensOIMergedCladereads.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(pathogenMergedCladereads, file = paste(outputDir,"/pathogens/",shortDate,"_pathogensOIMergedCladereads.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")

