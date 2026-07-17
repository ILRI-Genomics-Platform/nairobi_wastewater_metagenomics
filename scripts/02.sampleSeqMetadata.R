# Read wastewater sample metadata
#   - sample metadata
#   - sequence QC metadata
#   - classification metadata
# perform some cleaning
#   - remove duplicates
#   - remove specific runs: run08 & run09 - seems to have over-representation of some bacteterial pseudomonas species
# Plot bar plots of read counts
#   - All runs : raw read counts vs trimmed read counts
#   - Filtered runs :  raw read counts vs trimmed read counts

#Installing packages
setwd("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/")
lib="/home/gkibet/R/x86_64-pc-linux-gnu-library/4.3"
source("./scripts/functions.R")

requiredCRANPackages= c("dplyr", "ggplot2", "igraph", "readr", "openxlsx","tidytext",
                    "stringr", "xml2", "tidyverse", "magrittr", "webr",
                    "ggVennDiagram","devtools","svglite","janitor") #"ggvenn",
#  CRAN packages will be installed if they are not yet installed.
installNloadCRANpackages(requiredPackages = requiredCRANPackages, lib = lib)

# devtools::install_github("yanlinlin82/ggvenn@try-solve-count", force = TRUE)
# detach("package:ggvenn", unload = TRUE)
# library(ggvenn)

getwd() #Get working directory

#Loading data and wrangling 
#Loading metadata
dataDate="20230824"
sampleMetadata <- read.xlsx(paste0("./metadata/",dataDate,"_wastewater_samples_metadata.xlsx",sep=""), "sample_metadata", colNames = T, rowNames = F) %>%
  mutate(weekOfYear = strftime(as.Date(COLLECTION_DATE, format = "%d/%m/%Y"), format = "%V")) %>% 
  mutate(Year = strftime(as.Date(COLLECTION_DATE, format = "%d/%m/%Y"), format = "%y")) %>% 
  select(SAMPLE_NUMBER,weekOfYear,Year) %>% distinct()
siteMetadata <- read.csv(paste0("./metadata/",dataDate,"_wastewater_samples_metadata.csv",sep=""), header = T, row.names = NULL) %>%
  right_join(.,sampleMetadata)
siteCoordinates <- read.xlsx("./metadata/Nairobi_samplingSites_coordinates.xlsx", "sampling Sites", colNames = T, rowNames = F) %>%
  select(Estate,Latitude,Longitude)
seqMetadata <- read.xlsx(paste0("./metadata/",dataDate,"_wastewater_samples_metadata.xlsx",sep=""), "sequence_data", colNames = T, rowNames = F)
readQCMetadata <- read.csv(paste0("./plotdata/metrics/",dataDate,"_fastpQC_SummaryStatistics.csv",sep=""), header = T, sep = "\t")
readCounts <- readQCMetadata %>% select(- ends_with("bases"), raw.q30_rate) %>% 
  mutate(dropped.reads = raw.total_reads - trimmed.total_reads)

# Merging metadata
readMetadata <- left_join(seqMetadata,readCounts)
metadata <- right_join(siteMetadata,readMetadata, by = join_by(SAMPLE_NUMBER == Name)) %>%
  rename_all(., .funs = tolower) %>% rename(Name=sample_number, EstateOfOrigin=origin_estate,
                                            weekOfYear=weekofyear,Seqrun=seqrun,CountyOfOrigin=diagnosis_county)
#Finding duplicated sample Names
metadata %>% group_by(Name) %>% 
  filter(., Name %in% c(subset(.,duplicated(Name))$Name)) %>%
  arrange(Name) -> duplicateMetadata
#View(duplicateMetadata)
#colnames(metadata)

# Cleaning up metadata
# - drop run8 & run9 - abundancies are errortic
# - Remove duplicates retaining samples with highest reads
drop_rows <- c("") #c("run08","run09") #c("run08","run09") #
drop_names <- c("") # c("SPL032") 
drop_sampleids <- c("COVG00051_S2","COVG00090_S9","COVG0093_S14","COVG00156_S10") # c("SPL032") 
status = "all" #filtered|all
drop_EstateOfOrigin <- c("") #c("EASTLEIGH","KARIOBANGI") # c("") 
cleanMetadata <- metadata %>% .[!(.$Seqrun %in% drop_rows),] %>% .[!(.$Name %in% drop_names),] %>% 
  .[!(.$sampleid %in% drop_sampleids),] %>% arrange(desc(raw.total_reads)) %>% distinct(Name, .keep_all = T) 
cleanMetadata00 <- cleanMetadata %>% .[!(.$EstateOfOrigin %in% drop_EstateOfOrigin),] %>% 
  group_by(weekOfYear) %>% mutate(SampleFreqPerWeek=n_distinct(Name)) %>% ungroup() %>%
  group_by(EstateOfOrigin) %>% mutate(SampleFreqPerSite=n_distinct(Name)) %>% ungroup() %>%
  group_by(EstateOfOrigin,weekOfYear) %>% mutate(SampleFreqPerSitePerWeek=n_distinct(Name)) %>% ungroup() %>%
  mutate(year_Week= paste(year,weekOfYear, sep = "_"), .after = year) %>%
  group_by(year_Week) %>% mutate(weekNo = sprintf("week%02d",cur_group_id()), .after = year_Week) %>% ungroup()

shortDate <- gsub("-","",base::Sys.Date())
shortDate="20240604"
write.table(metadata, file = paste("./metadata/",shortDate,"_all_cleanSampleMetadata.csv", sep = ""),
            quote = F, row.names = FALSE, col.names= TRUE, sep = ",")
write.table(cleanMetadata00, file = paste("./metadata/",shortDate,"_cleanSampleMetadata.csv", sep = ""),
            quote = F, row.names = FALSE, col.names= TRUE, sep = ",")

plotData <- cleanMetadata00 %>%
  select(Name,EstateOfOrigin,Seqrun,trimmed.total_reads,dropped.reads,weekNo,CountyOfOrigin,
         SampleFreqPerWeek,SampleFreqPerSite,SampleFreqPerSitePerWeek) %>%
  pivot_longer(cols = ends_with("reads"), names_to = "readsCategory", values_to = "readsCount") %>%
  mutate_at("readsCategory", list(~str_replace(., "dropped.reads", "Trimmed"))) %>%
  mutate_at("readsCategory", list(~str_replace(., "trimmed.total_reads","Passed QC (=>20QScore)")))
plotData00 <- plotData %>% pivot_wider(names_from = readsCategory, values_from = readsCount) %>%
  select(-starts_with("Sample"))
write.table(plotData00, file = paste("./plotdata/metrics/",shortDate,"_fastpQC_SummaryStatistics_plotData.csv", sep = ""),
            quote = F, row.names = FALSE, col.names= TRUE, sep = ",")


# Plot Retained Plus trimmed reads Counts - faceted by EstateOfOrigin
ggplot(plotData, aes(y=readsCount, x=weekNo, fill = CountyOfOrigin, group = interaction(Name,readsCategory), 
                     colour = readsCategory)) +
  geom_col(position = position_stack( vjust = 1, reverse = T)) +
  theme(text = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right") +
  scale_y_continuous(labels = scales::comma) +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "CountyOfOrigin", x = "Sample ID", y = "Reads Count (Trimmed + dropped)") +
  scale_fill_viridis_d(option = "plasma") +
  facet_grid(cols = vars(EstateOfOrigin), scales = "free", space = "free")
ggsave(paste("./plots/metrics/",shortDate,"_",status,"_trimmed-DroppedReadCount-EstateofOrigin.png", sep = ""),
       width = 110, height = 20, units = "cm")

# Plot Retained Plus trimmed reads Counts - Stacked bar graph  - faceted by EstateOfOrigin
ggplot(plotData, aes(y=readsCount, x=weekNo, fill = CountyOfOrigin, group = interaction(Name,readsCategory),
                     colour = readsCategory)) +
  geom_col(position = position_fill( vjust = 1, reverse = T)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(text = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right") +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "CountyOfOrigin", x = "Sample ID", y = "Reads Count (Trimmed + dropped)") +
  scale_fill_viridis_d(option = "plasma") +
  facet_grid(. ~ EstateOfOrigin, scales = "free", space = "free")
ggsave(paste("./plots/metrics/",shortDate,"_",status,"_trimmed-DroppedReadCount-EstateofOrigin_stackedbar.png", sep = ""),
       width = 110, height = 20, units = "cm")

# Plot Retained Plus trimmed reads Counts - faceted by Weeks
ggplot(plotData, aes(y=readsCount, x=EstateOfOrigin, fill = CountyOfOrigin, group = interaction(Name,readsCategory), 
                     colour = readsCategory)) +
  geom_col(position = position_stack( vjust = 1, reverse = T)) +
  theme(text = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right") +
  scale_y_continuous(labels = scales::comma) +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "CountyOfOrigin", x = "Sample ID", y = "Reads Count (Trimmed + dropped)") +
  scale_fill_viridis_d(option = "plasma") +
  facet_grid(. ~ weekNo, scales = "free", space = "free")
ggsave(paste("./plots/metrics/",shortDate,"_",status,"_trimmed-DroppedReadCount-weekNo.png", sep = ""),
       width = 110, height = 30, units = "cm")

# Plot Retained Plus trimmed reads Counts - Stacked bar graph  - faceted by Weeks
ggplot(plotData, aes(y=readsCount, x=EstateOfOrigin, fill = CountyOfOrigin, group = interaction(Name,readsCategory),
                     colour = readsCategory)) +
  geom_col(position = position_fill( vjust = 1, reverse = T)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(text = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right") +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "CountyOfOrigin", x = "Sample ID", y = "Reads Count (Trimmed + dropped)") +
  scale_fill_viridis_d(option = "plasma") +
  facet_grid(. ~ weekNo, scales = "free", space = "free")
ggsave(paste("./plots/metrics/",shortDate,"_",status,"_trimmed-DroppedReadCount-weekNo_stackedbar.png", sep = ""),
       width = 110, height = 30, units = "cm")

# Plot Retained + trimmed read Counts facetted by runs
ggplot(plotData, aes(y=readsCount, x=Name, fill = CountyOfOrigin, group = interaction(Name,readsCategory),
                     colour = readsCategory)) +
  geom_bar(stat="identity", position = position_stack( vjust = 1, reverse = T)) +
  theme(text = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right") +
  scale_y_continuous(labels = scales::comma) +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "CountyOfOrigin", x = "Sample ID", y = "Reads Count (Passed QC 20QScore + dropped)") +
  scale_fill_viridis_d(option = "plasma") +
  facet_grid(. ~ Seqrun, scales = "free", space = "free") +
  geom_text(aes(label = EstateOfOrigin), vjust = 0.5, hjust = 0.0, angle = 90, size = 4, colour="black", check_overlap = TRUE) +
  geom_text(aes(label = readsCount), vjust = 0.5, hjust = 1, angle = 90, size = 4, colour="red", check_overlap = TRUE)
ggsave(paste("./plots/metrics/",shortDate,"_",status,"_trimmed-DroppedReadCount-RunID.png", sep = ""),
       width = 130, height = 30, units = "cm",limitsize = FALSE)

# Plot Retained + trimmed read Counts facetted by runs  - Stacked bar graph
ggplot(plotData, aes(y=readsCount, x=Name, fill = CountyOfOrigin, group = interaction(Name,readsCategory),
                     colour = readsCategory)) +
  geom_col(position = position_fill( vjust = 1, reverse = T)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))+
  theme(text = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right") +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "CountyOfOrigin", x = "Sample ID", y = "Reads Count (Trimmed + dropped)") +
  scale_fill_viridis_d(option = "plasma") +
  facet_grid(. ~ Seqrun, scales = "free", space = "free") +
  geom_text(aes(label = EstateOfOrigin), position = position_fill( vjust = 0.5, reverse = T), hjust = 1, angle = 90, size = 4, colour="black")
ggsave(paste("./plots/metrics/",shortDate,"_",status,"_trimmed-DroppedReadCount-RunID_stackedbar.png", sep = ""),
       width = 130, height = 30, units = "cm",limitsize = FALSE)

# Plot Raw read Counts
ggplot(cleanMetadata00, aes(y=raw.total_reads, x=Name, fill = CountyOfOrigin)) +
  geom_col(position = "stack") +
  theme(text = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right") +
  scale_y_continuous(labels = scales::comma) +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "CountyOfOrigin", x = "Sample ID", y = "Total Raw Reads Count") +
  scale_fill_viridis_d(option = "plasma") +
  facet_grid(. ~ EstateOfOrigin, scales = "free", space = "free")
ggsave(paste("./plots/metrics/",shortDate,"_",status,"_rawReadCount-EstateofOrigin.png", sep = ""),
       width = 110, height = 30, units = "cm")

ggplot(cleanMetadata00, aes(y=raw.total_reads, x=Name, fill = CountyOfOrigin)) +
  geom_col(position = "stack") +
  theme(text = element_text(face = "bold", size = 20), 
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.5),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right") +
  scale_y_continuous(labels = scales::comma) +
  guides(fill=guide_legend(ncol=1)) +
  labs(fill = "CountyOfOrigin", x = "Sample ID", y = "Total Raw Reads Count") +
  scale_fill_viridis_d(option = "plasma") +
  facet_grid(. ~ Seqrun, scales = "free", space = "free") +
  geom_text(aes(label = EstateOfOrigin), vjust = 0.5, hjust = 0.0, size = 7, angle = 90)
ggsave(paste("./plots/metrics/",shortDate,"_",status,"_rawReadCount-RunID.png", sep = ""),
       width = 110, height = 30, units = "cm")

