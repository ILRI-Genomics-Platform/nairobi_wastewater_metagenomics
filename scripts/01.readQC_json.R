#what this script does:
# 1. reads multiple fastp.json files
# 2. aggregate them into a list
# 3. Reads though the list, extracting 'summary' statistics and joins them all in one dataframe.
# 4. Calculates the mean for each summary statistic: total.reads, total.bases, ...

#Installation and setup
setwd("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/")
lib="/home/gkibet/R/x86_64-pc-linux-gnu-library/4.3"
source("./scripts/functions.R")

requiredCRANPackages= c("tidyverse", "rjson", "readr", "jsonlite", "tools",
                    "xml2", "magrittr","hrbrthemes","openxlsx")

#  CRAN packages will be installed if they are not yet installed.
installNloadCRANpackages(requiredPackages = requiredCRANPackages, lib = lib)

#Read from Sample_data.csv
runs=c("run00","run01","run02","run03","run04","run05","run00","run08","run09")
metadata <- read.xlsx("metadata/20230824_wastewater_samples_metadata.xlsx", colNames = T, 
                      sheet = "sequence_data") %>% filter(.,Seqrun %in% runs)

#View(metadata)
shortDate <- gsub("-","",base::Sys.Date())
#shortDate
#shortDate="20221116"

# Create input df for filter_Reports() - Must have two columns 'Name' and 'ReportFilePath'
# 'Name' column contains preferred names for the samples from which the report was generated
# 'ReportFilePath' column contains filepaths to the kreports matching the 'Name'
jsonFiles <- metadata %>% select(Name,sampleID,ReportFilePath) %>% 
  mutate(QCjsonFilePath = gsub('kreports/wastewater/centrifuge', 'QCreports',
                               gsub('.kreport.txt','.fastp.json', ReportFilePath, 
                                      ignore.case = FALSE, perl = FALSE, 
                                      fixed = FALSE, useBytes = FALSE)))
# path="/home/gkibet/bioinformatics/github/service/20240711_KEMRIWELLCOME_AA_NextSeqHT_SM_SANTHE6v2/results/fastp/"
# jsonFiles2 <- list.files(path = path, pattern = "*.json", full.names = T)
# jsonFiles3 <- jsonFiles2 %>% as.data.frame() %>% rename_at(1,~"QCjsonFilePath") %>%
#   mutate(sampleID = gsub( '.fastp', '', file_path_sans_ext(basename(QCjsonFilePath))))
# jsonFiles=jsonFiles3
#View(jsonFiles)

# Running the analysis
qcReportsList <- read_jsonQCReports(jsonFilesdf = jsonFiles,idColName="sampleID")
fastp_summary_df <- read_summaryQCreports(qcReportsList)

write.table(fastp_summary_df, file = paste("./plotdata/metrics/",shortDate,"_fastpQC_SummaryStatistics.csv", sep = ""),
            row.names = FALSE, col.names= TRUE, sep = "\t")
write.xlsx(fastp_summary_df, file = paste("./plotdata/metrics/",shortDate,"_fastpQC_SummaryStatistics.xlsx", sep = ""),
           rowNames = FALSE, colNames= TRUE, sep = "\t")

# readCounts <- fastp_summary_df %>% select(sampleID,raw.total_reads, trimmed.total_reads) %>% 
#   mutate(dropped.reads = raw.total_reads - trimmed.total_reads)

