#Installing packages 
rm(list = ls()) # List all objects in current environment and remove them
unlink(".Rhistory") # Detach and clear pre-existing R history
gc() # Free Unsused R memory
lib="/home/gkibet/R/x86_64-pc-linux-gnu-library/4.3"
# setwd("/home/gkibet/hpc_mnt/bioinformatics/github/metagenomics/data/visualization/")
# source("/home/gkibet/hpc_mnt/bioinformatics/github/metagenomics/data/visualization/scripts/functions.R")
setwd("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/")
source("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/scripts/functions.R")

# vector of packages to load
# CRAN
requiredCRANPackages= c("BiocManager","dplyr", "ggplot2", "readr", "openxlsx","stringr", "xml2", "tidyverse", "vroom",
                        "magrittr","ggVennDiagram","devtools","janitor","RColorBrewer","BBmisc","scales","Hmisc",
                        "igraph", "tidygraph", "graphlayouts", "ggraph","psych","RCy3","doParallel","foreach","parallel",
                        "data.table", "rjson", "readr", "jsonlite", "tools")
# BIOManager
requiredBIOCPackages=c("RCy3")
#  CRAN packages will be installed if they are not yet installed.
installNloadCRANpackages(requiredPackages = requiredCRANPackages, lib = lib)
# BiocManager and GitHub packages will be installed using the functions below if they are not yet installed. Uncomment to install.
# installNloadBIOCpackages(requiredPackages = requiredBIOCPackages, lib = lib)
# installNloadGitHubpackages(requiredPackages = requiredGITHUBPackages, lib = lib)
# getwd() #Get working directory

shortDate <- gsub("-","",base::Sys.Date())
shortDate="20241105"
shortDate="20250715"

# Metadata:
#Read from Sample_data.csv
runs=c("run00","run01","run02","run03","run04","run05","run00","run08","run09")
sampleSheet <- read.csv("../2023-06-12_run01-09_nextseq_metagen/samplesheet.csv",sep = ",",header = T) %>%
  select(sample,short_reads_1) %>% 
  mutate( Seqrun = str_extract(short_reads_1, "run\\d+"), sampleID = str_extract(short_reads_1, "COVG\\d+_S\\d+")) %>%
  select(-short_reads_1) %>% mutate(Seqrun = ifelse(Seqrun == "run35", "run00", Seqrun))
seq_data <- read.xlsx("metadata/20230824_wastewater_samples_metadata.xlsx", colNames = T, sheet = "sequence_data") %>% 
  filter(.,Seqrun %in% runs) %>% select(Name,sampleID,Seqrun)
seq_metadata <- merge(sampleSheet,seq_data, all.x = TRUE, all.y = TRUE) %>% select(Name,sample,sampleID,Seqrun)
sampleMetadata <- read.csv("./metadata/20230824_wastewater_samples_metadata.csv", sep = ",", header = T, row.names = NULL) %>%
  mutate(weekOfYear = strftime(as.Date(COLLECTION_DATE, format = "%d/%m/%Y"), format = "%V")) %>% 
  mutate(Year = strftime(as.Date(COLLECTION_DATE, format = "%d/%m/%Y"), format = "%y")) #%>% select(SAMPLE_NUMBER,weekOfYear,Year) %>% distinct()
rawMetadata <- sampleMetadata %>% merge(., seq_metadata, by.x = "SAMPLE_NUMBER", by.y="Name", all.y = TRUE)  %>%
  rename_all(., .funs = tolower) %>% 
  select(sample_number,origin_estate,diagnosis_county,collection_date,socioeconomic_category,weekofyear,year,sample,sampleid,seqrun) %>%
  rename(Name=sample_number, EstateOfOrigin=origin_estate,weekOfYear=weekofyear,Seqrun=seqrun,CountyOfOrigin=diagnosis_county) %>%
  mutate(year_Week= paste(year,weekOfYear, sep = "_"), .after = year) %>% select(-c("year","weekOfYear"))
  # group_by(year_Week) %>% mutate(weekNo = sprintf("week%02d",cur_group_id()), .after = year_Week) %>% ungroup()

# QC Reports:
# Fastp - read trimming and qc reports
# Create a list of Fastp output JSON files
file_paths <- list.files(path="results_arg/fastp",pattern = "COVG.*.fastp.json", full.names = T)
# Extract the names from the file paths
names <- sub(".*\\/([^\\/]+)\\.fastp\\.json$", "\\1", file_paths)
# Create a dataframe
jsonFiles <- data.frame(name = names, QCjsonFilePath = file_paths, stringsAsFactors = FALSE)
# Read the json files to a list of dataframes
fastpQCReportsList <- read_jsonFASTPReports(jsonFilesdf = jsonFiles,idColName="name")
# Merge the fastp dataframes from a list to a single df
fastp_df <- read_summaryFASTPQCreports(fastpQCReportsList)
# Write the fastp output
fastp_file = paste("./results_arg/r-analysis/data/",shortDate,"_QC_fastp.tsv", sep = "")
write.table(fastp_df, file = fastp_file, row.names = F, col.names = T, sep = "\t", quote = F)

# Hostile - host read removal and qc reports
# Create a list of Hostile output JSON files
hostileQCjson <- list.files(path="results_arg/hostile",pattern = "COVG.*.hostile.logs.json", full.names = T)
# Running the analysis
hostileReportsList <- read_jsonHostile(hostileQCjson)
# Merge the list of Hostile
hostile_df <- read_summaryHostileQCreports(hostileReportsList)
# Write the hostile output
hostile_file = paste("./results_arg/r-analysis/data/",shortDate,"_QC_hostile.tsv", sep = "")
write.table(hostile_df, file = hostile_file, row.names = F, col.names = T, sep = "\t", quote = F)

# Hostile + Fastp : merged QC reports
fastp_summary_df <- fastp_df %>% select(sampleID,raw.total_reads,raw.read1_mean_length,raw.read2_mean_length,
                                        trimmed.total_reads,trimmed.read1_mean_length,trimmed.read2_mean_length,trimmed.gc_content)
hostile_summary_df <- hostile_df %>% select(sampleID,reads_in,reads_out,reads_removed,reads_removed_proportion) %>%
  mutate(across(-1, ~ format(as.numeric(.), scientific = FALSE)))
qc_summary_df <- full_join(fastp_summary_df,hostile_summary_df, by = "sampleID")
# Merge all metadata
metadata <- qc_summary_df %>% select(sampleID,trimmed.total_reads) %>% 
  merge(., rawMetadata, by.x = "sampleID", by.y="sample", all.y = TRUE) %>%
  rename(sample_Name = sampleid)

# Write the hostile+fastp output
metadata_file = paste("./results_arg/r-analysis/data/",shortDate,"_sample_metadata.tsv", sep = "")
# use the format function to convert the numbers to character strings before writing them to a file - Avoids Exponentials
# qc_summary_df <- lapply(qc_summary_df, function(x) { if (is.numeric(x)) format(x, scientific = FALSE) else x})
write.table(metadata, file = metadata_file, row.names = F, col.names = T, sep = "\t", quote = F)
# metadata <- read.csv(metadata_file, sep = "\t", header = T)

# AMR Gene Families
# Define the patterns to remove 
patterns_to_remove <- c(" resistance", " \\(5-nitroimidazole\\)", " \\(Glycopeptid resistance\\)")
resFinderAMRClasses <- data.frame(Line = read_lines(file = "https://bitbucket.org/genomicepidemiology/resfinder_db/raw/e0525f203e43bcd6dc91348e561a2cb09fd4014b/notes.txt"),
             stringsAsFactors = FALSE ) %>%
  mutate(AMRClass = ifelse(grepl("^#", Line), sub(" ","_",sub(" *resistance","",sub(":","",sub("^# *", "", Line)))), NA)) %>%
  fill(AMRClass, .direction = "down") %>%
  filter(!grepl("^#", Line)) %>%
  separate(Line, into = c("Gene.Name", "DrugClass", "alternateGeneName"), sep = ":") %>%
  mutate(DrugClass = str_replace_all(DrugClass, setNames(rep("", length(patterns_to_remove)), patterns_to_remove)))
fileName="./results_arg/r-analysis/data/resFinder_AMRGeneClasses_DrugClass_Unique.tsv"
write.table(resFinderAMRClasses, file = fileName, row.names = F, col.names = T,
            sep = "\t", quote = F, na = "NA")
  
fileName="./rgi/args/resFinderCARD_AMRGeneFamily_DrugClass_ARGs_Unique.csv"
AMRGeneFamilies <- read.csv(fileName,sep = ",",header = T,na.strings = c("", "NA")) %>%
  select(AMR.Gene.Family,AMR.Gene.Family.Name,resFinder.Gene.Name) %>% 
  filter(!is.na(resFinder.Gene.Name)) %>%
  separate_rows(resFinder.Gene.Name, sep = ";")

# ARG datasets:
# Create a list of ResFinder output JSON files
resfinderQCjson <- list.files(path="results_arg/resfinder",pattern = "COVG.*.resfinder.json", full.names = T)
# Running the analysisq
resFinderReportsList <- read_jsonResfinder(resfinderQCjson)
# Merge the list of resFinder AMR classes/Phenotypes
AMRDrugList <- filter_ResfinderReports(resFinderReportsList)
# Create a dataframe of unique ref_id and name
seqRegionsAMR_DF <- AMRDrugList$seqRegionsDF %>% 
  select(ref_id,name,ref_acc,ref_seq_length) %>% 
  distinct() %>% as.data.frame()
# Write the hostile output
AMR_file = paste("./results_arg/r-analysis/data/",shortDate,"_SeqRegions_AMR.tsv", sep = "")
phenotype_file = paste("./results_arg/r-analysis/data/",shortDate,"_Phenotypes_Drugs.tsv", sep = "")
write.table(AMRDrugList$seqRegionsDF, file = AMR_file, row.names = F, col.names = T, sep = "\t", quote = F)
write.table(AMRDrugList$phenotypesDF, file = phenotype_file, row.names = F, col.names = T, sep = "\t", quote = F)
# AMRDrugList <- list()
# AMRDrugList$seqRegionsDF <- read.csv(AMR_file, sep = "\t", header = T)
# AMRDrugList$phenotypesDF <- read.csv(phenotype_file, sep = "\t", header = T)

# Read Counts from ARGs:
# Create a list of ResFinder output JSON files
resfinderfraggz <- list.files(
  path="results_arg/resfinder",pattern = "kma_.*.frag.gz", full.names = T, recursive = T)
# Merge the list of resfinder frag.gz files
outputList <- read_resFinderFragGz(resfinderfraggz)
AMRAbundance_df <- outputList$AMRAbundance_df
# Write the hostile output
AMRAbundance_file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRAbundance_df.tsv", sep = "")
write.table(AMRAbundance_df, file = AMRAbundance_file, row.names = F, col.names = T, sep = "\t", quote = F)
# AMRAbundance_df <- read.csv(AMRAbundance_file, sep = "\t", header = T)

# Merge AMRAbundance_df (read fragments data) and seqRegionsAMR_DF (AMRClass occurence data) and filter reads <100 mappingScore
seqRegionsDF <- merge(AMRAbundance_df, seqRegionsAMR_DF, by.x = "templateID", by.y = "ref_id", all.y = TRUE) %>%
  select(everything(), templateID, readID, name) 
seqRegionsDF00 <- seqRegionsDF %>% filter(mappingScore >= 100)
# seqRegionsDF$name %>% is.na() %>% table() %>% as.data.frame()
# name_iddf <- seqRegionsDF %>% filter(is.na(name)) %>% select(templateID,name) %>% distinct() %>% as.data.frame()
# name_id_df <- seqRegionsDF %>% filter(!is.na(name)) %>% select(templateID,name) %>% distinct() %>% as.data.frame()
# seqRegionsAMR_DF1 <- seqRegionsAMR_DF %>% arrange(name)
# setdiff(seqRegionsAMR_DF1$name,name_id_df$name)
# setdiff(name_id_df$name, seqRegionsAMR_DF1$name)
# setdiff(seqRegionsAMR_DF1$ref_id,name_id_df$templateID)
# setdiff(name_id_df$templateID, seqRegionsAMR_DF1$ref_id)

# Merge forward and reverse read records into one
AMRAbundance_AMRClass_df <- seqRegionsDF %>% 
  arrange(readID) %>%
  mutate(readID_base = sub("/[12]$", "", readID)) %>%
  group_by(readID_base, sampleID, databaseID, drugClass, 
           templateID, name, ref_acc, ref_seq_length) %>%
  do(merge_rows(.)) %>%
  ungroup() %>%
  select(-readID_base)
# Write the hostile output
AMRAbundance_AMRClass_file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRAbundance_edited_df.tsv", sep = "")
# AMRAbundance_AMRClass_file = paste("./results_arg/r-analysis/data/20241105_AMRAbundance_edited_df0.tsv", sep = "")
write.table(AMRAbundance_AMRClass_df, file = AMRAbundance_AMRClass_file, row.names = F, col.names = T, sep = "\t", quote = F)
# AMRAbundance_AMRClass_df <- read.csv(AMRAbundance_AMRClass_file, sep = "\t", header = T)

AMRAbundance_AMRClass_df00 <- merge(
  AMRAbundance_AMRClass_df,AMRGeneFamilies,
  by.x = "name", by.y = "resFinder.Gene.Name", all.x = TRUE)
# outdfAMR <- AMRAbundance_AMRClass_df00 %>% 
#   select(sampleID,name,amrClass,templateID,drugClass,AMR.Gene.Family) %>% 
#   distinct() %>% filter(is.na(AMR.Gene.Family)) #%>% filter(is.na(name))

# Count abundances per sampleID+databaseID, sampleID+databaseID+drugClass, sampleID+databaseID+drugClass+templateID
AMRAbundance_AMRClass_NoNA <- AMRAbundance_AMRClass_df00 %>% ungroup() %>%
  filter(!is.na(amrClass), databaseID=="resfinder") %>% 
  group_by(sampleID,databaseID) %>%
  mutate(readAbundancePerDatabase = n_distinct(readID)) %>% ungroup() %>%
  group_by(sampleID,drugClass) %>%
  mutate(readAbundancePerDrugClass = n_distinct(readID)) %>% ungroup() %>%
  group_by(sampleID,amrClass) %>%
  mutate(readAbundancePerAMRClass = n_distinct(readID)) %>% ungroup() %>%
  group_by(sampleID,AMR.Gene.Family) %>%
  mutate(readAbundancePerAMRGeneFamily = n_distinct(readID)) %>% ungroup() %>%
  group_by(sampleID,templateID) %>%
  mutate(readAbundancePerTemplateID = n_distinct(readID)) %>% ungroup()

# Write the hostile output
ARGAbundance_Nmetadata_file = paste("./results_arg/r-analysis/data/",shortDate,"_ARGAbundance_And_metadata.tsv", sep = "")
write.table(AMRAbundance_AMRClass_NoNA, file = ARGAbundance_Nmetadata_file, row.names = F, col.names = T, sep = "\t", quote = F)
# AMRAbundance_AMRClass_NoNA <- read.csv(ARGAbundance_Nmetadata_file, sep = "\t", header = T)


# Merge AMRAbundance_AMRClass_NoNA and metadata
AMRAbundance_Nmetadata_df <- metadata %>% #filter(sampleID == Name) %>% 
  select(Name,sampleID,trimmed.total_reads,EstateOfOrigin,collection_date,socioeconomic_category,year_Week,Seqrun) %>%
  merge(., AMRAbundance_AMRClass_NoNA, by.x = "Name", by.y="sampleID", all.x = TRUE) %>%
  mutate(relativeAbundancePerDatabase = readAbundancePerDatabase/(trimmed.total_reads/2)) %>%
  mutate(relativeAbundancePerDrugClass = readAbundancePerDrugClass/(trimmed.total_reads/2)) %>%
  mutate(relativeAbundancePerAMRClass = readAbundancePerAMRClass/(trimmed.total_reads/2)) %>%
  mutate(relativeAbundancePerAMRGeneFamily = readAbundancePerAMRGeneFamily/(trimmed.total_reads/2)) %>%
  mutate(relativeAbundancePerTemplateID = readAbundancePerTemplateID/(trimmed.total_reads/2))  %>%
  mutate(RPK_TemplateID = readAbundancePerTemplateID/(ref_seq_length / 1000) ) %>%
  mutate(RPKM_TemplateID = ((RPK_TemplateID/(trimmed.total_reads/2)) * 1000000 )) %>%
  group_by(sampleID,amrClass) %>%
  mutate(RPK_AMRClass = sum(RPK_TemplateID)) %>%
  mutate(RPKM_AMRClass = ((RPK_AMRClass/(trimmed.total_reads/2)) * 1000000 )) %>%
  group_by(sampleID,drugClass) %>%
  mutate(RPK_drugClass = sum(RPK_TemplateID)) %>%
  mutate(RPKM_drugClass = ((RPK_drugClass/(trimmed.total_reads/2)) * 1000000 )) %>%
  group_by(sampleID) %>%
  mutate(TPM_TemplateID = 1e6 * RPK_TemplateID / sum(RPK_TemplateID, na.rm = TRUE)) %>%
  mutate(TPM_AMRClass = 1e6 * RPK_AMRClass / sum(RPK_TemplateID, na.rm = TRUE)) %>%
  mutate(TPM_drugClass = 1e6 * RPK_drugClass / sum(RPK_TemplateID, na.rm = TRUE)) %>%
  ungroup()
#%>% select(-c("read"))
# Write the hostile output
AMRAbundance_Nmetadata_file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRAbundance_Nmetadata.tsv", sep = "")
# AMRAbundance_Nmetadata_file = paste("./results_arg/r-analysis/data/20241105_AMRAbundance_Nmetadata0.tsv", sep = "")
write.table(AMRAbundance_Nmetadata_df, file = AMRAbundance_Nmetadata_file, row.names = F, col.names = T, sep = "\t", quote = F)
# AMRAbundance_Nmetadata_df <- read.csv(AMRAbundance_Nmetadata_file, sep = "\t", header = T)
AMRGeneFamily_metadata_df <- AMRAbundance_Nmetadata_df %>% 
  select(name,drugClass,amrClass,AMR.Gene.Family,AMR.Gene.Family.Name) %>% distinct()
AMRGeneFamily_metadata_file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRGeneFamily_metadata.tsv", sep = "")
write.table(AMRGeneFamily_metadata_df, file = AMRGeneFamily_metadata_file, row.names = F, col.names = T, sep = "\t", quote = F)

# Now let us explore the seq_regions metrics:
colnames(AMRDrugList$seqRegionsDF)
histIdentity <- hist(AMRDrugList$seqRegionsDF$identity)
histCoverage <- hist(AMRDrugList$seqRegionsDF$coverage)
# Write plots to file:
# In pdf format
# Open a PDF device
pdf(paste("./results_arg/r-analysis/plots/",shortDate,"_Identity_Coverage_histograms.pdf", sep = ""))
# Set up the plotting area to have 2 rows and 1 column
par(mfrow = c(2, 1))
# Plot the histograms
plot(histIdentity)
plot(histCoverage)
# Close the PDF device
dev.off()

# In pdf format
# Open a JPEG device
jpeg(paste("./results_arg/r-analysis/plots/",shortDate,"_Identity_Coverage_histograms.jpeg", sep = ""))
# Set up the plotting area to have 2 rows and 1 column
par(mfrow = c(2, 1))
# Plot the histograms
plot(histIdentity)
plot(histCoverage)
# Close the PDF device
dev.off()

vplotData=seqRegionsDF %>% mutate(mappingCategory = ifelse(mappingScore <= 100, "1-100", "Over 100"))
ggplot(vplotData, aes(x=sampleID, y=mappingScore, fill=sampleID)) + 
  geom_violin() +
  geom_jitter(aes(color = mappingCategory)) +
  theme(text = element_text(face = "bold", size = 30),
        axis.text.y = element_text(angle = 90),
        legend.position="none",
        plot.title = element_text("mappingScore")) +
  ylab("mappingScore") +
  scale_fill_viridis_d(option = "plasma") +
  scale_color_manual(values = c("1-100" = "orange", "Over 100" = "gray"))+
  facet_wrap(~ sampleID, nrow = 4, scales = "free_x")
fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_sampleID_mappingScore_violinPlot.png", sep = "")
ggsave(fileName, width = 200, height = 100, units = "cm",limitsize = FALSE)

# Create the plot
ggplot(vplotData, aes(x = sampleID, y = mappingScore, fill = sampleID)) + 
  geom_violin() +
  geom_boxplot(width = 0.2) +
  theme(text = element_text(face = "bold", size = 30),
        axis.text.y = element_text(angle = 90),
        legend.position = "none",
        plot.title = element_text("mappingScore")) +
  ylab("mappingScore") +
  scale_fill_viridis_d(option = "plasma") +
  scale_color_manual(values = c("1-100" = "orange", "Over 100" = "gray")) +
  facet_wrap(~ sampleID, nrow = 4, scales = "free_x")
fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_sampleID_mappingScore_boxplot.png", sep = "")
ggsave(fileName, width = 200, height = 100, units = "cm",limitsize = FALSE)
