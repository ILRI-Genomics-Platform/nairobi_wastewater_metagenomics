#Installing packages
setwd("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/")
rm(list = ls()) # List all objects in current environment and remove them
unlink(".Rhistory") # Detach and clear pre-existing R history
gc() # Free Unsused R memory
lib="/home/gkibet/R/x86_64-pc-linux-gnu-library/4.3"
source("/home/gkibet/bioinformatics/github/metagenomics/data/visualization/scripts/functions.R")

# vector of packages to load
# CRAN
requiredCRANPackages= c("BiocManager","dplyr", "ggplot2", "readr", "openxlsx","stringr", "xml2", "tidyverse","magrittr",
                        "ggVennDiagram","devtools","janitor","RColorBrewer","BBmisc","scales","Hmisc","igraph", "tidygraph",
                        "graphlayouts", "ggraph","psych","RCy3","data.table", "rjson", "readr", "jsonlite", "tools",
                        "vegan", "ggrepel","dendextend","pairwiseAdonis")
# BIOManager
requiredBIOCPackages=c("RCy3","impute")
# GitHub
requiredGITHUBPackages=c("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#  CRAN packages will be installed if they are not yet installed.
installNloadCRANpackages(requiredPackages = requiredCRANPackages, lib = lib)
# BiocManager packages will be installed using the functions below if they are not yet installed. Uncomment to install.
installNloadBIOCpackages(requiredPackages = requiredBIOCPackages, lib = lib)
# GitHub packages will be installed using the functions below if they are not yet installed. Uncomment to install.
# installNloadGitHubpackages(requiredPackages = requiredGITHUBPackages, lib = lib)
# getwd() #Get working directory

shortDate <- gsub("-","",base::Sys.Date())
shortDate="20241105"

# Read the AMR Class - Drug Class Abundance datasets
drop_EstateOfOrigin = c("EASTLEIGH","KARIOBANGI")
drop_Seqrun = c("run09","run08")
drop_sampleIDs = c("COVG090b","COVG093b","COVG156b")
AMRAbundance_Nmetadata_file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRAbundance_Nmetadata.tsv", sep = "")
AMRAbundance_Nmetadata_df0 <- read.csv(AMRAbundance_Nmetadata_file, sep = "\t", header = T)
AMRAbundance_Nmetadata_df <- AMRAbundance_Nmetadata_df0 %>% 
  filter(!EstateOfOrigin %in% drop_EstateOfOrigin) %>% #,"DANDORA"
  filter(!Seqrun %in% drop_Seqrun)

# Now let us explore the abundance metrics:
colnames(AMRAbundance_Nmetadata_df)
# Generate relative abundance matrices of drugClass(phenotyes) and amrClass(genes)

# Get unique drug classes
drug_classes <- unique(AMRAbundance_Nmetadata_df$drugClass)

# Initialize an empty list to store matrices
matrix_list <- list()
# Loop through each drug class and create the matrix
for (drug_Class in drug_classes) {
  cat("\tProcessing data from",drug_Class,"Drug class\n")
  # Function to create the matrix for a given drugClass
  create_matrix <- function(drug_Class, data) {
    data %>%
      filter(drugClass == !!drug_Class) %>%
      select(sampleID, amrClass, relativeAbundancePerAMRClass) %>%
      distinct() %>%
      pivot_wider(names_from = amrClass, values_from = relativeAbundancePerAMRClass) %>%
      as.data.frame() %>%
      column_to_rownames(var = "sampleID")
  }
  
  matrix <- create_matrix(drug_Class, AMRAbundance_Nmetadata_df)
  matrix_list[[drug_Class]] <- matrix
}

# Now matrix_list contains a matrix for each drugClass
print(matrix_list)

# Initialize an empty data frame to store the metadata
metadata_df <- data.frame(
  MatrixName = character(),
  NumObservations = integer(),
  NumVariables = integer(),
  Variables = character(),
  stringsAsFactors = FALSE
)

# Loop through the list of matrices and extract metadata
for (name in names(matrix_list)) {
  # name="tetracycline"
  cat("\tProcessing data matrix from",name,"Drug class\n")
  # Function to extract metadata from a matrix
  extract_metadata <- function(matrix, name) {
    num_obs <- nrow(matrix)
    num_vars <- ncol(matrix)
    
    # Calculate colMeans(is.na(<column>)) for each column
    na_means <- colMeans(is.na(matrix)) %>% sort()
    # Create a string with column names and their NA means
    vars_with_na <- sapply(names(na_means), function(v) paste(v, "=", round(na_means[v], 3), sep = ""))
    vars_str <- paste(vars_with_na, collapse = ";")
    
    data.frame(
      MatrixName = name,
      NumObservations = num_obs,
      NumVariables = num_vars,
      Variables = vars_str,
      stringsAsFactors = FALSE
    )
  }
  
  matrix <- matrix_list[[name]]
  metadata <- extract_metadata(matrix, name)
  metadata_df <- bind_rows(metadata_df, metadata)
}

# View the resulting metadata data frame
print(metadata_df)
write.table(metadata_df, row.names = F, col.names = T, sep = "\t", quote = F,
            file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRDrugClassMatrices_Metadata.tsv", sep = ""))

# Metadata
inMetadata00 <- AMRAbundance_Nmetadata_df %>% select(sampleID,EstateOfOrigin,collection_date,socioeconomic_category,year_Week,Seqrun) %>%
  distinct() %>% arrange(year_Week) %>% group_by(year_Week) %>% 
  mutate(weekNo = sprintf("week%02d",cur_group_id()), .after = year_Week) %>% #ungroup()
  mutate(socioeconomic_category = if_else(EstateOfOrigin == "DANDORA", "mixed", socioeconomic_category)) %>% ungroup()
  # mutate(socioeconomic_category = if_else(socioeconomic_category == "middle_income", "low_income", socioeconomic_category)) %>% ungroup()
# 20241105 - DANDORA considered mixed(socialeconomic_category)
# 20241105b - DANDORA included in low_income
# 20241105b1 - DANDORA included in low_income, looking into Seqruns
# 20241105c - DANDORA dropped
# 20241105d - grouped by seq_run - looking for run08 & run09
# # Exclude run08 qnd run09
# 20241120 - DANDORA considered mixed(socialeconomic_category)
# 20241120b - DANDORA included in low_income
# 20241120c - DANDORA dropped
# 20241120d - DANDORA dropped + All middle_income as high_income
# 20241120e - DANDORA dropped + All middle_income as low_income

shortDate="20241120"
# relativeAbundancePerDrugClass_df <- relativeAbundancePerDrugClass_mtx %>% 
#   rownames_to_column(var = "sampleID")
# write.table(merge(inMetadata,relativeAbundancePerDrugClass_df), row.names = F, col.names = T, sep = "\t", quote = F,
#             file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRData_Nmetadata.tsv", sep = ""))

# Create a log file and initialise it with an empty string.
logFileName=paste("./results_arg/r-analysis/OrdinationPERMANOVA_Analysis_",format(Sys.time(), "%Y%m%d_%H%M%S"),".log", sep="")
writeLines(c(""),logFileName)
# Loop through the list of matrices and perform NMDS and PCoA based on distance matrices:
for (AMRDrugClass in names(matrix_list)) {
  # Redirect R output to the log file connection. Appends console output.
  sink(file = logFileName, append=TRUE)
  # sink()
  # Input data test
  # AMRDrugClass = "colistin"
  # AMRDrugClass = "beta-lactam"
  # AMRDrugClass = "tetracycline"
  inPutData = matrix_list[[AMRDrugClass]]
  inPutDateName = paste("relativeAbundancePerAMRClass_",AMRDrugClass,"Mtx",sep = "")
  inMetadata <- inMetadata00 %>% filter(sampleID %in% row.names(inPutData))
  
  # Remove columns with more than a certain percentage of NAs
  threshold <- 0.6  # Adjust this threshold as needed
  inPutData_filtered <- inPutData[, colMeans(is.na(inPutData)) < threshold]
  AMRGeneCount = ncol(inPutData_filtered)
  inPutData_filteredMerged <- merge(rownames_to_column(inPutData_filtered, var = "sampleID"),inMetadata) %>% 
    group_by(EstateOfOrigin,socioeconomic_category) %>% 
    mutate(samplesPerEstateOfOrigin = n()) %>% ungroup()
  fileName=paste("./results_arg/r-analysis/data/",shortDate,"_",inPutDateName,
                 "_Filtered_Genes_below",gsub("0\\.","Point",threshold),"NAmeans_",AMRGeneCount,".tsv", sep = "")
  write.table(inPutData_filteredMerged, row.names = F, col.names = T, sep = "\t", quote = F,
              file = fileName)
  
  # Bray-Curtis dissimilarity test
  # Bray-Curtis dissimilarity test matrix
  # dissimilarity_algorthmn <- "BrayCurtis_dissimilarity"
  bray_curtis_matrix <- vegdist(inPutData, method = "bray", na.rm = TRUE)
  write.table(as.matrix(bray_curtis_matrix), 
              file = paste("./results_arg/r-analysis/data/",shortDate,"_bray_curtis_",inPutDateName,".tsv", sep = ""),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # Aitchison distance (also known as the Aitchison distance or CLR distance)
  # Aitchison distance test matrix
  # dissimilarity_algorthmn <- "RobustAitchison" #"RobustAitchison"|"Aitchison"
  robust_aitchison_matrix <- vegdist(inPutData, method = "robust.aitchison", na.rm = TRUE)
  write.table(as.matrix(robust_aitchison_matrix), 
              file = paste("./results_arg/r-analysis/data/",shortDate,"_robust_aitchison_",inPutDateName,".tsv", sep = ""),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # Assign dissimilarity matrix
  distanceMatrixList <- list(BrayCurtis_dissimilarity = bray_curtis_matrix,
                             RobustAitchison = robust_aitchison_matrix)
  for (distanceMatrix in names(distanceMatrixList)) {
    # distanceMatrix = "BrayCurtis_dissimilarity"
    # distanceMatrix = "RobustAitchison"
    dissimilarity_matrix00 <- distanceMatrixList[[distanceMatrix]]
    dissimilarity_algorthmn <- distanceMatrix 
    cat("\tProcessing",dissimilarity_algorthmn,"data matrix from",AMRDrugClass,"Drug class\n")
    
    # Perform imputation using the impute package
    dissimilarity_matrix <- as.dist(impute.knn(as.matrix(dissimilarity_matrix00))$data)
    
    # Perform hierarchical clustering 
    hc <- hclust(dissimilarity_matrix, method = "ward.D2")
    dend <- as.dendrogram(hc, center = T)
    plot(dend, main = "Hierarchical Clustering Dendrogram", ylab = "Height",
         xlab= paste(dissimilarity_algorthmn,"matrix",sep = " "))
    
    # Ordination
    # Perform PCoA (cmdscale())
    cat("\t\tPCOA analysis on",dissimilarity_algorthmn,AMRDrugClass,"Drug class matrix\n")
    dimentionality_reduction <- "PCoA"
    pcoa_results <- cmdscale(dissimilarity_matrix, eig = TRUE, k = 2)
    # Convert PCoA results to a dataframe
    pcoa_df0 <- as.data.frame(pcoa_results$points) %>% rename_all(~str_replace_all(., "V", "PCoA"))
    pcoa_df <- merge(rownames_to_column(pcoa_df0, var = "sampleID"),inMetadata)
    pcoA_Plot_base <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2)) + 
      geom_point(aes(color = EstateOfOrigin, shape = socioeconomic_category), size = 3) + 
      labs(title = "PCoA Plot", x = "PCoA Axis 1", y = "PCoA Axis 2") + 
      theme(text = element_text(face = "plain", size = 20),
            axis.text.y = element_text(angle = 0), legend.position="right")
    pcoA_Plot_base_ellipse <- pcoA_Plot_base +
      stat_ellipse(aes(colour = socioeconomic_category), linewidth = 1) + # linewidth not size, since ggplot 3.4.0
      coord_fixed(ratio = 1, clip = "off")
    # linetype = socioeconomic_category,
    
    # Calculate the species scores that represent the contribution of each species to the principal coordinates.
    # envfit function from the vegan package fits environmental variables (species) to the ordination
    pcoa_points <- pcoa_df0 %>% select(PCoA1, PCoA2) %>% as.matrix() %>% na.omit(.)
    
    # Impute remaining NAs
    inPutData_imputed <- setNames(
      as.data.frame(lapply(inPutData_filtered, function(x) {
        if (is.numeric(x)) {
          x[is.na(x)] <- mean(x, na.rm = TRUE)
        }
        return(x)
      })),
      names(inPutData_filtered)
    )
    # inPutData_imputed = inPutData_filtered
    
    envfit_results <- envfit(pcoa_results, inPutData_imputed, perm = 999, na.rm = TRUE)
    # Extract the scores from the envfit object and convert them to a data frame
    amrClass_scores <- as.data.frame(scores(envfit_results, display = "vectors")) %>%
      rename_all(~str_replace_all(., "Dim", "PCoA"))
    amrClass_scores$amrClass <- rownames(amrClass_scores)
    # Calculate the length of vectors to determine importance 
    amrClass_scores$length <- sqrt(amrClass_scores$PCoA1^2 + amrClass_scores$PCoA2^2)
    amrClass_scores <- amrClass_scores %>% 
      mutate(labelAngle = atan2(PCoA2, PCoA1) * 180 / pi, label_x = PCoA1 * 1.1, label_y = PCoA2 * 1.1)
    # Filter for the most important vectors (e.g., top 10) 
    top_n <- 100 # Adjust this number as needed 
    important_vectors <- amrClass_scores %>% arrange(desc(length)) %>% head(top_n)
    
    # Calculate the scaling factor: PCoA1 and PCoA2 axes are much larger than the lengths of vectors
    # Causes vectors to appear very short and congested at the center.
    # Need to rescale the vectors to better match the scale of your PCoA axes
    division_fastor <- max(abs(c(pcoa_df$PCoA1, pcoa_df$PCoA2)))
    scaling_factor <- max(abs(c(important_vectors$PCoA1, important_vectors$PCoA2))) / division_fastor
    important_vectors$PCoA1_scaled <- important_vectors$PCoA1 / scaling_factor
    important_vectors$PCoA2_scaled <- important_vectors$PCoA2 / scaling_factor
    # Add the vectors to the plot: Use geom_segment to add the vectors and geom_text to label them
    pcoA_Plot <- pcoA_Plot_base_ellipse +
      geom_segment(data = important_vectors, aes(x = 0, y = 0, xend = PCoA1_scaled, yend = PCoA2_scaled), 
                   arrow = arrow(length = unit(0.2, "cm")), color = "grey50") +
      geom_label_repel(data = important_vectors, aes(x = PCoA1_scaled, y = PCoA2_scaled, label = amrClass, angle = labelAngle), 
                       fill = "white", color = "black", box.padding = unit(0.3, "lines"), 
                       point.padding = unit(0.3, "lines"), segment.color = "black", 
                       size = 2) +
      annotate("text", x = Inf, y = Inf, label = paste("Bi-plots rescaled by division by", scaling_factor),
               hjust = 1, vjust = 1, size = 4, color = "grey50") + 
      guides(color = guide_legend(override.aes = list(shape = 19, size = 3)), # Override shapes in legend
             shape = guide_legend(override.aes = list(size = 3))) # Ensure shape sizes are consistent
    
    # Display the plot
    pcoA_Plot
    
    # Lets Save the plot now
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",dissimilarity_algorthmn,
                   "_important_vectors_Top",top_n,"_PCoA_Plot.png", sep = "")
    dirPath <- dirname(fileName)
    # Create the directory if it does not exist
    if (!dir.exists(dirPath)) { dir.create(dirPath, recursive = TRUE) }
    
    ggsave(fileName, plot = pcoA_Plot, width = 35, height = 25, units = "cm",limitsize = FALSE)
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",dissimilarity_algorthmn,
                   "_important_vectors_Top",top_n,"_PCoA_Plot.pdf", sep = "")
    ggsave(fileName, plot = pcoA_Plot, width = 35, height = 25, units = "cm",limitsize = FALSE,
           device = cairo_pdf)
    
    # Perform NMDS
    cat("\t\tNMDS analysis on",dissimilarity_algorthmn,AMRDrugClass,"Drug class matrix\n")
    dimentionality_reduction <- "NMDS"
    nmds_results <- metaMDS(dissimilarity_matrix, k = 2, trymax = 100)
    # Convert PCoA results to a dataframe
    nmds_df0 <- as.data.frame(nmds_results$points) %>% rename_all(~str_replace_all(., "MDS", "NMDS"))
    nmds_df <- merge(rownames_to_column(nmds_df0, var = "sampleID"),inMetadata)
    nmds_Plot_base <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) + 
      geom_point(aes(color = EstateOfOrigin, shape = socioeconomic_category), size = 3) + 
      labs(title = "NMDS Plot", x = "NMDS Axis 1", y = "NMDS Axis 2") + 
      theme(text = element_text(face = "plain", size = 20),
            axis.text.y = element_text(angle = 0),
            legend.position="right") 
    nmds_Plot_base_ellipse <- nmds_Plot_base +
      stat_ellipse(aes(colour = socioeconomic_category), linewidth = 1) + # linewidth not size, since ggplot 3.4.0
      coord_fixed(ratio = 1, clip = "off")
    # linetype = socioeconomic_category, 
    # Calculate the species scores that represent the contribution of each species to the principal coordinates.
    # envfit function from the vegan package fits environmental variables (species) to the ordination
    nmds_points <- nmds_df0 %>% select(NMDS1, NMDS2) %>% as.matrix() %>% na.omit(.)
    
    # Calculate environmental vectors
    envfit_results <- envfit(nmds_results, inPutData_imputed, perm = 999, na.rm = TRUE)
    # Extract the scores from the envfit object and convert them to a data frame
    amrClass_scores <- as.data.frame(scores(envfit_results, display = "vectors"))
    amrClass_scores$amrClass <- rownames(amrClass_scores)
    # Calculate the length of vectors to determine importance 
    amrClass_scores$length <- sqrt(amrClass_scores$NMDS1^2 + amrClass_scores$NMDS2^2)
    amrClass_scores <- amrClass_scores %>% 
      mutate(angle = atan2(NMDS1, NMDS2) * 180 / pi, label_x = NMDS1 * 1.1, label_y = NMDS2 * 1.1)
    # Filter for the most important vectors (e.g., top 10) 
    top_n <- 100 # Adjust this number as needed 
    important_vectors <- amrClass_scores %>% arrange(desc(length)) %>% head(top_n)
    
    # Calculate the scaling factor: PCoA1 and PCoA2 axes are much larger than the lengths of vectors
    # Causes vectors to appear very short and congested at the center.
    # Need to rescale the vectors to better match the scale of your PCoA axes
    division_fastor <- max(abs(c(nmds_df$NMDS1, nmds_df$NMDS2)))
    scaling_factor <- max(abs(c(important_vectors$NMDS1, important_vectors$NMDS2))) / division_fastor
    important_vectors$NMDS1_scaled <- important_vectors$NMDS1 / scaling_factor
    important_vectors$NMDS2_scaled <- important_vectors$NMDS2 / scaling_factor
    
    # Add the vectors to the plot: Use geom_segment to add the vectors and geom_text to label them
    nmds_Plot <- nmds_Plot_base_ellipse +
      geom_segment(data = important_vectors, aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled), 
                   arrow = arrow(length = unit(0.2, "cm")), color = "grey50") +
      geom_label_repel(data = important_vectors, aes(x = NMDS1_scaled, y = NMDS2_scaled, label = amrClass), 
                       fill = "white", color = "black", box.padding = unit(0.3, "lines"), 
                       point.padding = unit(0.3, "lines"), segment.color = "black", 
                       size = 2, max.overlaps = 30) +
      annotate("text", x = Inf, y = Inf, label = paste("Bi-plots rescaled by division by", scaling_factor),
               hjust = 1, vjust = 1, size = 4, color = "grey50")
    
    # Display the plot
    nmds_Plot
    
    # Lets Save the plot now
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",
                   dissimilarity_algorthmn,"_important_vectors_Top",top_n,"_NMDS_Plot.png", sep = "")
    ggsave(fileName, plot = nmds_Plot, width = 35, height = 25, units = "cm",limitsize = FALSE)
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",
                   dissimilarity_algorthmn,"_important_vectors_Top",top_n,"_NMDS_Plot.pdf", sep = "")
    ggsave(fileName, plot = nmds_Plot, width = 35, height = 25, units = "cm",limitsize = FALSE)
    
    
    # PERMANOVA (Permutational Multivariate Analysis of Variance)
    # Performing PERMANOVA permutes/shaffles the lables* many times (~999) to create a distribution of differences that
    # counld happen by chance. Then determine if the observed differences between groups are greater than what would be 
    # expected by random chance. It is a powerful way to test for significant differences in community composition among groups.
    # Helps explain if different groups (e.g., habitats, treatments, or time points) have significantly different 
    # compositions of species or other multivariate data.
    # How do we do PERMANOVA: use "adonis()" function of the "vegan" package in R.
    
    # Model: Include socioeconomic_category, EstateOfOrigin, weekNo, and Seqrun as predictors.
    # by = "margin": Tests the marginal effect of each variable by permuting each variable independently of others
    # Goal: assess the unique contribution of each variable to the dissimilarity.
    cat("\tPerforming PERMANOVA on",dissimilarity_algorthmn,"data matrix from",AMRDrugClass,"Drug class\n")
    permanova_results <- adonis2(dissimilarity_matrix ~ socioeconomic_category + EstateOfOrigin + weekNo + Seqrun, 
                                 data = inMetadata, permutations = 999, by = "margin")
    print(permanova_results)
    
    # Save the results to a text file
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",
                   dissimilarity_algorthmn,"_PERMANOVA_results.txt", sep = "")
    writeLines(capture.output((permanova_results)), con = fileName, sep = "\n", useBytes = TRUE)
    
    # Model: Similar to (1) but excludes EstateOfOrigin as a predictor.
    # by = "margin": Tests the marginal effect of each variable without considering EstateOfOrigin.
    # Goal: exclude EstateOfOrigin. socioeconomic_category is a grouping of and EstateOfOrigin*, they rule out each other.
    permanova_results0 <- adonis2(dissimilarity_matrix ~ socioeconomic_category + weekNo + Seqrun, 
                                  data = inMetadata, permutations = 999, by = "margin")
    print(permanova_results0)
    
    # Save the results to a text file
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",
                   dissimilarity_algorthmn,"_PERMANOVA_margin.txt", sep = "")
    writeLines(capture.output((permanova_results0)), con = fileName, sep = "\n", useBytes = TRUE)
    
    # Model: Includes socioeconomic_category, EstateOfOrigin, weekNo, and Seqrun as predictors.
    # by = "terms": Tests the sequential (type I) sums of squares. Tests the addition of each term one by one in the order specified.
    # assess the effect of each variable accounting the variables listed before it in the model.
    permanova_results1 <- adonis2(dissimilarity_matrix ~ socioeconomic_category + EstateOfOrigin + weekNo + Seqrun, 
                                  data = inMetadata, permutations = 999, by = "terms")
    print(permanova_results1)
    
    # Save the results to a text file
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",
                   dissimilarity_algorthmn,"_PERMANOVATerms_results.txt", sep = "")
    writeLines(capture.output((permanova_results1)), con = fileName, sep = "\n", useBytes = TRUE)
    
    
    # Model: Includes only socioeconomic_category as a predictor.
    # by (default is "terms"): Tests the overall effect of socioeconomic_category alone.
    permanova_results2 <- adonis2(dissimilarity_matrix ~ socioeconomic_category, 
                                  data = inMetadata, permutations = 999)
    print(permanova_results2)
    
    # Save the results to a text file
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",
                   dissimilarity_algorthmn,"_PERMANOVA_SCWk_results.txt", sep = "")
    writeLines(capture.output((permanova_results2)), con = fileName, sep = "\n", useBytes = TRUE)
    
    # Model: Includes only socioeconomic_category as a predictor.
    # strata = inMetadata$weekNo: Stratifies permutations by weekNo, permutations are constrained within weekNo strata.
    # useful as there is a hierarchical structure in the data ~ accounts for the variability within each week.
    permanova_results3 <- adonis2(dissimilarity_matrix ~ socioeconomic_category, 
                                  data = inMetadata, permutations = 999,
                                  strata = inMetadata$weekNo)
    print(permanova_results3)
    
    # Save the results to a text file
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",
                   dissimilarity_algorthmn,"_PERMANOVA_strataweekNo_results.txt", sep = "")
    writeLines(capture.output((permanova_results3)), con = fileName, sep = "\n", useBytes = TRUE)
    
    # Pairwise adonis : pairwiseAdonis package
    # Assuming 'socioeconomic_category' is a factor indicating the group (low_income, middle_income, high_income)
    pairwiseAdonis_results <- pairwise.adonis(dissimilarity_matrix, inMetadata$socioeconomic_category)
    print(pairwiseAdonis_results)
    
    # Save the results to a text file
    fileName=paste("./results_arg/r-analysis/plots/",AMRDrugClass,"/",shortDate,"_",inPutDateName,"_",
                   dissimilarity_algorthmn,"_pairwiseAdonis_results.txt", sep = "")
    writeLines(capture.output((pairwiseAdonis_results)), con = fileName, sep = "\n", useBytes = TRUE)
    
    # Stop diverting output to the log file
    sink() 
  }
}

