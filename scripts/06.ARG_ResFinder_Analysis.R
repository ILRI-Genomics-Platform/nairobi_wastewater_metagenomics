#Installing packages
setwd("/home/gkibet/hpc_mnt/bioinformatics/github/metagenomics/data/visualization/")
rm(list = ls()) # List all objects in current environment and remove them
unlink(".Rhistory") # Detach and clear pre-existing R history
gc() # Free Unsused R memory
lib="/home/gkibet/R/x86_64-pc-linux-gnu-library/4.4"
source("/home/gkibet/hpc_mnt/bioinformatics/github/metagenomics/data/visualization/scripts/functions.R")

# vector of packages to load
# CRAN
requiredCRANPackages= c(
  "BiocManager", "dplyr", "ggplot2", "readr", "openxlsx", "stringr", "xml2", "tidyverse",
  "magrittr", "ggVennDiagram", "devtools", "janitor", "RColorBrewer", "BBmisc", "scales", 
  "Hmisc", "igraph", "tidygraph", "graphlayouts", "ggraph", "psych", "RCy3", "data.table",
  "rjson", "readr", "jsonlite", "tools", "vegan", "ggrepel", "dendextend", "pairwiseAdonis",
  "ggvenn", "pheatmap", "webr", "ggdendro", "gtable", "svglite", "grid", "proxy", 
  "RColorBrewer", "gridExtra", "patchwork", "ggplotify")

# BIOManager
requiredBIOCPackages=c("RCy3")
# GitHub
requiredGITHUBPackages=c("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#  CRAN packages will be installed if they are not yet installed.
installNloadCRANpackages(requiredPackages = requiredCRANPackages, lib = lib)
# BiocManager and GitHub packages will be installed using the functions below if they are not yet installed. Uncomment to install.
# installNloadBIOCpackages(requiredPackages = requiredBIOCPackages, lib = lib)
# installNloadGitHubpackages(requiredPackages = requiredGITHUBPackages, lib = lib)
# getwd() #Get working directory

shortDate <- gsub("-","",base::Sys.Date())
shortDate="20241105"
# shortDate="20260715"

# Read the AMR Class - Drug Class Abundance datasets
# Pre-filtered data:
AMRAbundance_Nmetadata_file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRAbundance_Nmetadata.tsv", sep = "")
AMRAbundance_Nmetadata_df0 <- read.csv(AMRAbundance_Nmetadata_file, sep = "\t", header = T)
AMRAbundance_Nmetadata_df <- AMRAbundance_Nmetadata_df0 %>% 
  filter(!EstateOfOrigin %in% c("EASTLEIGH","KARIOBANGI")) %>% #,"DANDORA"
  filter(!Seqrun %in% c("run09","run08")) %>%
  mutate(
    socioeconomic_category = case_when(
      EstateOfOrigin == "DANDORA" ~ "central_WWTP",
      TRUE ~ socioeconomic_category))

# Unfiltered data:
AMRAbundance_Nmetadata_file = paste("./results_arg/r-analysis/data/20241105_AMRAbundance_Nmetadata0.tsv", sep = "")
AMRAbundance_resfinder <- read.csv(AMRAbundance_Nmetadata_file, sep = "\t", header = T) %>%
  group_by(name) %>% mutate(ARO.Term.Freq = n()) %>% ungroup() %>%
  filter(!EstateOfOrigin %in% c("EASTLEIGH","KARIOBANGI"))  %>%
  mutate(
    socioeconomic_category = case_when(
      EstateOfOrigin == "DANDORA" ~ "central_WWTP",
      TRUE ~ socioeconomic_category)) %>% 
  as.data.frame()

AMRAbundance_resfinder_uniqueARGs <- AMRAbundance_resfinder %>% 
  select("databaseID", "drugClass", "amrClass") %>% 
  distinct()

# Now let us explore the abundance metrics:
colnames(AMRAbundance_Nmetadata_df)
colnames(AMRAbundance_resfinder)
length(unique(AMRAbundance_Nmetadata_df$amrClass))
length(unique(AMRAbundance_Nmetadata_df$drugClass))
length(unique(AMRAbundance_Nmetadata_df$name))
length(unique(AMRAbundance_Nmetadata_df$Seqrun))
length(unique(AMRAbundance_Nmetadata_df$Name))
length(unique(AMRAbundance_resfinder$amrClass))
length(unique(AMRAbundance_resfinder$drugClass))
length(unique(AMRAbundance_resfinder$name))
length(unique(AMRAbundance_resfinder$Seqrun))
length(unique(AMRAbundance_resfinder$Name))

# Generate relative abundance matrices of drugClass(phenotyes) and amrClass(genes)
# Extract relativeAbundancePerAMRClass
cat("relativeAbundancePerAMRClass_mtx abundance matrices...\n")
relativeAbundancePerAMRClass_mtx <- AMRAbundance_Nmetadata_df %>% 
  select(sampleID,amrClass,relativeAbundancePerAMRClass) %>%
  distinct() %>% as.data.frame() %>% 
  pivot_wider(names_from = amrClass, 
              values_from = relativeAbundancePerAMRClass) %>%
  as.data.frame(lapply(., as.numeric)) %>% 
  column_to_rownames(var = "sampleID")

# Extract RPKMPerAMRClass
cat("RPKMPerAMRClass_mtx abundance matrices...\n")
RPKMPerAMRClass_mtx <- AMRAbundance_Nmetadata_df %>% 
  select(sampleID,amrClass,RPKM_AMRClass) %>%
  distinct() %>% as.data.frame() %>% 
  pivot_wider(names_from = amrClass, 
              values_from = RPKM_AMRClass) %>%
  as.data.frame(lapply(., as.numeric)) %>% 
  column_to_rownames(var = "sampleID")

# Extract TPMPerAMRClass
cat("TPMPerAMRClass_mtx abundance matrices...\n")
TPMPerAMRClass_mtx <- AMRAbundance_Nmetadata_df %>% 
  select(sampleID,amrClass,TPM_AMRClass) %>%
  distinct() %>% as.data.frame() %>% 
  pivot_wider(names_from = amrClass, 
              values_from = TPM_AMRClass) %>%
  as.data.frame(lapply(., as.numeric)) %>% 
  column_to_rownames(var = "sampleID")

# Extract readAbundancePerDrugClass
relativeAbundancePerDrugClass_mtx <- AMRAbundance_Nmetadata_df %>% 
  select(sampleID,drugClass,relativeAbundancePerDrugClass) %>%
  distinct() %>% as.data.frame() %>% 
  pivot_wider(names_from = drugClass, 
              values_from = relativeAbundancePerDrugClass) %>%
  as.data.frame(lapply(., as.numeric)) %>% 
  column_to_rownames(var = "sampleID")

# Extract RPKMPerDrugClass
cat("RPKMPerDrugClass_mtx abundance matrices...\n")
RPKMPerDrugClass_mtx <- AMRAbundance_Nmetadata_df %>% 
  select(sampleID,drugClass,RPKM_drugClass) %>%
  distinct() %>% as.data.frame() %>% 
  pivot_wider(names_from = drugClass, 
              values_from = RPKM_drugClass) %>%
  as.data.frame(lapply(., as.numeric)) %>% 
  column_to_rownames(var = "sampleID")

# Extract TPMPerDrugClass
cat("TPMPerDrugClass_mtx abundance matrices...\n")
TPMPerDrugClass_mtx <- AMRAbundance_Nmetadata_df %>% 
  select(sampleID,drugClass,TPM_drugClass) %>%
  distinct() %>% as.data.frame() %>% 
  pivot_wider(names_from = drugClass, 
              values_from = TPM_drugClass) %>%
  as.data.frame(lapply(., as.numeric)) %>% 
  column_to_rownames(var = "sampleID")

# Extract relativeAbundancePerAMRGeneFamily
relativeAbundancePerAMRGeneFamily_mtx <- AMRAbundance_Nmetadata_df %>% 
  select(sampleID,AMR.Gene.Family,relativeAbundancePerAMRGeneFamily) %>%
  distinct() %>% as.data.frame() %>% 
  pivot_wider(names_from = AMR.Gene.Family, 
              values_from = relativeAbundancePerAMRGeneFamily) %>%
  as.data.frame(lapply(., as.numeric)) %>% 
  column_to_rownames(var = "sampleID")

# Metadata
inMetadata <- AMRAbundance_Nmetadata_df %>% 
  select(sampleID,EstateOfOrigin,collection_date,socioeconomic_category,year_Week,Seqrun) %>%
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

# shortDate="20241120"
shortDate="20260716"

relativeAbundancePerDrugClass_df <- relativeAbundancePerDrugClass_mtx %>% 
  rownames_to_column(var = "sampleID")
write.table(merge(inMetadata,relativeAbundancePerDrugClass_df), row.names = F, col.names = T, sep = "\t", quote = F,
            file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRData_Nmetadata.tsv", sep = ""))
write.table(AMRAbundance_resfinder_uniqueARGs, row.names = F, col.names = T, sep = "\t", quote = F,
            file = paste("./results_arg/r-analysis/data/",shortDate,"_AMRAbundance_resfinder_uniqueARGs.tsv", sep = ""))

# ####################################################################
cat("Saving abundace matrices...\n")
resfinderdatabasePath="./results_arg/r-analysis/data/"

# write AMRAbundance_Nmetadata_df to file
relativeAbundance_file <- paste(
  resfinderdatabasePath,shortDate,
  "_AMRData_relativeAbundance_Nmetadata.tsv", sep = "")
cat("Saving",relativeAbundance_file,"data...\n")
write.table(merge(inMetadata,AMRAbundance_Nmetadata_df), 
            row.names = F, col.names = T, sep = "\t", quote = F,
            file = relativeAbundance_file)
# AMRAbundance_Nmetadata_df <- read.csv(relativeAbundance_file, sep="\t", header=T)
# str(relativeAbundance_df)

# write relativeAbundancePerAMRClass_mtx to file
relativeAbundancePerAMRClass_df <- relativeAbundancePerAMRClass_mtx %>% 
  rownames_to_column(var = "sampleID")
relativeAbundancePerAMRClass_file <- paste(
  resfinderdatabasePath,shortDate,
  "_AMRData_relativeAbundancePerAMRClass_Nmetadata.tsv", sep = "")
cat("Saving",relativeAbundancePerAMRClass_file,"matrix...\n")
write.table(merge(inMetadata,relativeAbundancePerAMRClass_df), 
            row.names = F, col.names = T, sep = "\t", quote = F,
            file = relativeAbundancePerAMRClass_file)
# relativeAbundancePerAMRClass_mtx <- read.csv(
#   relativeAbundancePerAMRClass_file, sep = "\t", header = T) %>% 
#   column_to_rownames(var = "sampleID")
# str(relativeAbundancePerAMRClass_df)

# write RPKMPerAMRClass_mtx to file
RPKMPerAMRClass_df <- RPKMPerAMRClass_mtx %>% 
  rownames_to_column(var = "sampleID")
RPKMPerAMRClass_file <- paste(
  resfinderdatabasePath,shortDate,
  "_AMRData_RPKMPerAMRClass_Nmetadata.tsv", sep = "")
cat("Saving",RPKMPerAMRClass_file,"matrix...\n")
write.table(merge(inMetadata,RPKMPerAMRClass_df), 
            row.names = F, col.names = T, sep = "\t", quote = F,
            file = RPKMPerAMRClass_file)
# RPKMPerAMRClass_mtx <- read.csv(
#   RPKMPerAMRClass_file, sep = "\t", header = T) %>% 
#   column_to_rownames(var = "sampleID")
# str(RPKMPerAMRClass_df)

# write TPMPerAMRClass_mtx to file
TPMPerAMRClass_df <- TPMPerAMRClass_mtx %>% 
  rownames_to_column(var = "sampleID")
TPMPerAMRClass_file <- paste(
  resfinderdatabasePath,shortDate,
  "_AMRData_TPMPerAMRClass_Nmetadata.tsv", sep = "")
cat("Saving",TPMPerAMRClass_file,"matrix...\n")
write.table(merge(inMetadata,TPMPerAMRClass_df), 
            row.names = F, col.names = T, sep = "\t", quote = F,
            file = TPMPerAMRClass_file)
# TPMPerAMRClass_mtx <- read.csv(
#   TPMPerAMRClass_file, sep = "\t", header = T) %>% 
#   column_to_rownames(var = "sampleID")
# str(TPMPerAMRClass_df)

# write relativeAbundancePerDrugClass_mtx to file
relativeAbundancePerDrugClass_df <- relativeAbundancePerDrugClass_mtx %>% 
  rownames_to_column(var = "sampleID")
relativeAbundancePerDrugClass_file <- paste(
  resfinderdatabasePath,shortDate,
  "_AMRData_relativeAbundancePerDrugClass_Nmetadata.tsv", sep = "")
cat("Saving",relativeAbundancePerDrugClass_file,"matrix...\n")
write.table(merge(inMetadata,relativeAbundancePerDrugClass_df), 
            row.names = F, col.names = T, sep = "\t", quote = F,
            file = relativeAbundancePerDrugClass_file)
# relativeAbundancePerDrugClass_mtx <- read.csv(
#   relativeAbundancePerDrugClass_file, sep = "\t", header = T) %>% 
#   column_to_rownames(var = "sampleID")
# str(relativeAbundancePerDrugClass_df)

# write RPKMPerDrugClass_mtx to file
RPKMPerDrugClass_df <- RPKMPerDrugClass_mtx %>% 
  rownames_to_column(var = "sampleID")
RPKMPerDrugClass_file <- paste(
  resfinderdatabasePath,shortDate,
  "_AMRData_RPKMPerDrugClass_Nmetadata.tsv", sep = "")
cat("Saving",RPKMPerDrugClass_file,"matrix...\n")
write.table(merge(inMetadata,RPKMPerDrugClass_df), 
            row.names = F, col.names = T, sep = "\t", quote = F,
            file = RPKMPerDrugClass_file)
# RPKMPerDrugClass_mtx <- read.csv(
#   RPKMPerDrugClass_file, sep = "\t", header = T) %>% 
#   column_to_rownames(var = "sampleID")
# str(RPKMPerDrugClass_df)

# write TPMPerDrugClass_mtx to file
TPMPerDrugClass_df <- TPMPerDrugClass_mtx %>% 
  rownames_to_column(var = "sampleID")
TPMPerDrugClass_file <- paste(
  resfinderdatabasePath,shortDate,
  "_AMRData_TPMPerDrugClass_Nmetadata.tsv", sep = "")
cat("Saving",TPMPerDrugClass_file,"matrix...\n")
write.table(merge(inMetadata,TPMPerDrugClass_df), 
            row.names = F, col.names = T, sep = "\t", quote = F,
            file = TPMPerDrugClass_file)
# TPMPerDrugClass_mtx <- read.csv(
#   TPMPerDrugClass_file, sep = "\t", header = T) %>% 
#   column_to_rownames(var = "sampleID")
# str(TPMPerDrugClass_df)

# write relativeAbundancePerAMRGeneFamily_mtx to file
relativeAbundancePerAMRGeneFamily_df <- relativeAbundancePerAMRGeneFamily_mtx %>% 
  rownames_to_column(var = "sampleID")
relativeAbundancePerAMRGeneFamily_file <- paste(
  resfinderdatabasePath,shortDate,
  "_AMRData_relativeAbundancePerAMRGeneFamily_Nmetadata.tsv", sep = "")
cat("Saving",relativeAbundancePerAMRGeneFamily_file,"matrix...\n")
write.table(merge(inMetadata,relativeAbundancePerAMRGeneFamily_df), 
            row.names = F, col.names = T, sep = "\t", quote = F,
            file = relativeAbundancePerAMRGeneFamily_file)
# relativeAbundancePerAMRGeneFamily_mtx <- read.csv(
#   relativeAbundancePerAMRGeneFamily_df, sep = "\t", header = T) %>% 
#   column_to_rownames(var = "sampleID")
# str(relativeAbundancePerAMRGeneFamily_df)

# write inMetadata to file
inMetadata_file <- paste(
  resfinderdatabasePath,shortDate,"_AMRData_Nmetadata.tsv", sep = "")
cat("Saving",inMetadata_file,"matrix...\n")
write.table(inMetadata, row.names = F, col.names = T, sep = "\t", 
            quote = F, file = inMetadata_file)
# inMetadata <- read.csv(inMetadata_file, sep = "\t", header = T)
# str(inMetadata)

cat("Done saving relative abundance matrix data...\n")

# ####################################################################

## ============================================================================
## Reproduce: AMR Class heatmap (Z-scores, clustered, site-level) + 
##            3-way Venn diagram of gene overlap across income categories
## ============================================================================
## Assumes AMRAbundance_Nmetadata_df (filtered, site-collapsed heatmap input)
## and AMRAbundance_resfinder (unfiltered, used only for the Venn gene sets)
## already exist in your environment, as built in your earlier code.
## ============================================================================
# ---------------------------------------------------------------------------
# 0. Parameters -- adjust these to match your actual selection criteria
# ---------------------------------------------------------------------------
abundance_col <- "relativeAbundancePerAMRClass"   # relativeAbundancePerAMRClass or "RPKM_AMRClass"
n_top_classes <- 50                               # set NA to keep all amrClass
agg_fun       <- mean                             # how replicate samples per site are combined

# ---------------------------------------------------------------------------
# 1. Collapse sample-level data to SITE level (Name = estate),
#    since the heatmap's columns are sites, not individual sampleIDs
# ---------------------------------------------------------------------------
site_meta <- AMRAbundance_Nmetadata_df %>%
  distinct(EstateOfOrigin, socioeconomic_category)

site_amr <- AMRAbundance_Nmetadata_df %>%
  distinct(sampleID, EstateOfOrigin, amrClass, drugClass, !!sym(abundance_col)) %>%
  group_by(EstateOfOrigin, amrClass, drugClass) %>%
  summarise(value = agg_fun(.data[[abundance_col]], na.rm = TRUE), .groups = "drop")

# ---------------------------------------------------------------------------
# 2. Optional: keep only the top-N amrClass by overall mean abundance
#    (drop this filter, or replace with your real selection rule, e.g. a
#    prevalence cutoff, if you have a different original criterion)
# ---------------------------------------------------------------------------
if (!is.na(n_top_classes)) {
  top_classes <- site_amr %>%
    group_by(amrClass) %>%
    summarise(overall_mean = mean(value, na.rm = TRUE)) %>%
    arrange(desc(overall_mean)) %>%
    slice_head(n = n_top_classes) %>%
    pull(amrClass)
  
  site_amr <- site_amr %>% filter(amrClass %in% top_classes)
}

# ---------------------------------------------------------------------------
# 3. Build the site x amrClass matrix
# ---------------------------------------------------------------------------
amr_wide <- site_amr %>%
  select(EstateOfOrigin, amrClass, value) %>%
  pivot_wider(names_from = EstateOfOrigin, values_from = value) %>%
  column_to_rownames("amrClass") %>%
  as.matrix()

# write amr_wide to file
amr_wide_df <- as.data.frame(amr_wide) %>%
  rownames_to_column(var = "amrClass")
amr_wide_file <- paste(
  resfinderdatabasePath,shortDate,"_Top",n_top_classes,"_ARGs_Heatmap.tsv", sep = "")
cat("Saving",amr_wide_file,"matrix...\n")
write.table(amr_wide_df, row.names = F, col.names = T, sep = "\t", 
            quote = F, file = amr_wide_file)
# amr_wide_df <- read.csv(inMetadata_file, sep = "\t", header = T)
# str(amr_wide_df)

# ---------------------------------------------------------------------------
# 4. Row-wise Z-scores (each amrClass scaled across sites) -- this is what's
#    both plotted as the heatmap color AND printed as the cell text
# ---------------------------------------------------------------------------
z_mtx <- t(scale(t(amr_wide)))

# ---------------------------------------------------------------------------
# 5. Row annotation: Drug Class (one color-coded strip on the left)
# ---------------------------------------------------------------------------
row_anno <- site_amr %>%
  distinct(amrClass, drugClass) %>%
  column_to_rownames("amrClass")
row_anno <- row_anno[rownames(z_mtx), , drop = FALSE]   # keep row order matched

# ---------------------------------------------------------------------------
# 6. Column annotation: Socio-economic category (top strip)
# ---------------------------------------------------------------------------
col_anno <- site_meta %>%
  filter(EstateOfOrigin %in% colnames(z_mtx)) %>%
  distinct(EstateOfOrigin, .keep_all = TRUE) %>%
  column_to_rownames("EstateOfOrigin")
col_anno <- col_anno[colnames(z_mtx), , drop = FALSE]

# ---------------------------------------------------------------------------
# 7. Colors
# ---------------------------------------------------------------------------
drug_classes <- sort(unique(row_anno$drugClass))
drug_class_colors <- setNames(
  c("#F4A6C6", "#F08080", "#00CFCF", "#B8860B", "#E066B0",
    "#2E8B22", "#87CEFA", "#008B8B", "#C8A2C8", "#B19CD9", "#FFA500")[seq_along(drug_classes)],
  drug_classes
)

socio_colors <- c(
  "high_income"   = "#F4A6D6",
  "low_income"    = "#7CB342",
  "middle_income" = "#5C6BC0",
  "central_WWTP"  = "#9E9E9E"
)
socio_colors <- socio_colors[unique(col_anno$socioeconomic_category)]

annotation_colors <- list(
  drugClass = drug_class_colors,
  socioeconomic_category = socio_colors
)

# Diverging blue-white-red palette, centered at 0, matching the ~ -1 to 2 range
z_min <- floor(min(z_mtx, na.rm = TRUE))
z_max <- ceiling(max(z_mtx, na.rm = TRUE))
breaks <- c(seq(z_min, 0, length.out = 50), seq(0.01, z_max, length.out = 50))
heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks) - 1)

# ---------------------------------------------------------------------------
# 8. Build the heatmap
# ---------------------------------------------------------------------------
heatmap_plot <- pheatmap(
  z_mtx,
  color             = heat_colors,
  breaks            = breaks,
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  clustering_method = "complete",
  display_numbers   = matrix(sprintf("%.2f", z_mtx), nrow = nrow(z_mtx)),
  number_color      = "black",
  fontsize_number   = 6,
  na_col            = "grey80",
  annotation_row    = row_anno,
  annotation_col    = col_anno,
  annotation_colors = annotation_colors,
  border_color      = "white",
  fontsize_row      = 8,
  fontsize_col      = 9,
  angle_col         = 90,
  main              = "Antibiotic Resistance Genes Relative Abundance (Z-Scores)",
  silent            = TRUE
)

# grid.newpage()
# grid.draw(heatmap_plot$gtable)

# ---------------------------------------------------------------------------
# 9. Panel B: 3-way Venn of AMR gene (amrClass) presence across income levels
#    Built from the UNFILTERED resfinder table, excluding the central WWTP,
#    which is why totals sum to 207 (matches your unique(amrClass) count)
# ---------------------------------------------------------------------------
gene_sets <- AMRAbundance_resfinder %>% 
  mutate(socioeconomic_category = case_when(
    EstateOfOrigin == "DANDORA" ~ "low_income",
    TRUE ~ socioeconomic_category)) %>%
  filter(socioeconomic_category %in% c("low_income", "middle_income", "high_income")) %>%
  distinct(socioeconomic_category, amrClass) %>%
  group_by(socioeconomic_category) %>%
  summarise(genes = list(unique(amrClass)), .groups = "drop")

gene_list <- setNames(gene_sets$genes, gene_sets$socioeconomic_category)
# re-order to match figure legend layout: low, middle, high
gene_list <- gene_list[c("low_income", "middle_income", "high_income")]

venn_plot <- ggvenn(
  gene_list,
  fill_color    = c("#7CB342", "#5C6BC0", "#F4A6D6"),
  stroke_size   = 0.5,
  set_name_size = 6,
  text_size     = 4, 
  show_percentage = TRUE
)

print(venn_plot)

# ---------------------------------------------------------------------------
# 9b. Gene-level set differences, to list under the Venn diagram
# ---------------------------------------------------------------------------
low    <- gene_list[["low_income"]]
mid    <- gene_list[["middle_income"]]
high   <- gene_list[["high_income"]]
 
low_only        <- setdiff(low, union(mid, high))                    # unique to low-income
low_mid_only    <- setdiff(intersect(low, mid), high)                 # low & middle, not high
low_high_only   <- setdiff(intersect(low, high), mid)                 # low & high, not middle
 
# sanity check against the Venn counts (6 / 6 / 2 in your figure)
cat("low-only:", length(low_only),
    " | low & middle only:", length(low_mid_only),
    " | low & high only:", length(low_high_only), "\n")
 
gene_list_text <- paste0(
  "Unique to low-income (n=", length(low_only), "): ",
  paste(sort(low_only), collapse = ", "), "\n\n",
  "Shared: low-income & middle-income only (n=", length(low_mid_only), "): ",
  paste(sort(low_mid_only), collapse = ", "), "\n\n",
  "Shared: low-income & high-income only (n=", length(low_high_only), "): ",
  paste(sort(low_high_only), collapse = ", ")
)
 
# a small ggplot "panel" purely for displaying wrapped text under the Venn
gene_list_panel <- ggplot() +
  theme_void() +
  annotate(
    "text",
    x = 0, y = 1,
    label = stringr::str_wrap(gene_list_text, width = 55),
    hjust = 0, vjust = 1,
    size  = 3.2,
    lineheight = 1.1
  ) +
  xlim(0, 1) + ylim(0, 1)
 
# ---------------------------------------------------------------------------
# 10. Combine both panels into ONE figure, matching the source PDF layout:
#     Panel A = Z-score heatmap (left, wide)
#     Panel B = gene-overlap Venn diagram (top-right, smaller)
# ---------------------------------------------------------------------------

# pheatmap returns a grid grob, not a ggplot -- convert it so patchwork can use it
heatmap_gg <- as.ggplot(heatmap_plot)

# --- Add a black box around each panel individually ---
box_theme <- theme(
  plot.background = element_rect(colour = "black", fill = NA, linewidth = 1),
  plot.margin     = margin(8, 8, 8, 8)
)

heatmap_gg_boxed <- heatmap_gg + box_theme
venn_gene_panel  <- (venn_plot / gene_list_panel) +
  plot_layout(heights = c(2.2, 1.3)) # Venn gets more vertical space than the text
venn_gene_panel_boxed <- wrap_elements(full = venn_gene_panel) + box_theme

combined_plot <- (heatmap_gg_boxed | venn_gene_panel_boxed) +
  plot_layout(widths = c(2.4, 1)) +      # heatmap gets ~2.4x the width of panel B
  plot_annotation(tag_levels = "A")       # auto-labels panels "A" and "B"

print(combined_plot)

# Save at a size large enough for the dense heatmap text to stay legible
fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",abundance_col,"_Top",
               "",n_top_classes,"_AMR_heatmap_venn_combined.pdf", sep = "")
ggsave(
  filename = fileName,
  plot     = combined_plot,
  width    = 14,
  height   = 10,
  units    = "in",
  device = cairo_pdf
)

fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",abundance_col,"_Top",
               "",n_top_classes,"_AMR_heatmap_venn_combined.png", sep = "")
ggsave(
  filename = fileName,
  plot     = combined_plot,
  width    = 14,
  height   = 10,
  units    = "in",
  dpi      = 600
)

## ============================================================================
## Notes / things to double-check against your original analysis:
##  - abundance_col: switch to "RPKM_AMRClass" if that's what was actually used
##  - n_top_classes / selection rule: your figure shows 49 gene rows; the exact
##    original cutoff (top-N by abundance vs. a prevalence threshold) isn't
##    recoverable from the image alone -- adjust the filter in step 2 as needed
##  - agg_fun: mean() across replicate samples per site is a reasonable default;
##    swap to sum() if your original workflow re-derived relative abundance
##    from summed raw counts instead of averaging per-sample proportions
## ============================================================================

# # Input data test
# inPutData=relativeAbundancePerAMRGeneFamily_mtx #relativeAbundancePerDrugClass_mtx|relativeAbundancePerAMRClass_mtx|relativeAbundancePerAMRGeneFamily_mtx
# inPutDateName="relativeAbundancePerAMRGeneFamily_mtx" #"relativeAbundancePerDrugClass_mtx"|"relativeAbundancePerAMRClass_mtx"|"relativeAbundancePerAMRGeneFamily_mtx"
# drugClass="" #tetracycline|nitroimidazole|beta-lactam|macrolide|phenicol|sulphonamide|trimethoprim|quinolone|rifampicin|aminoglycoside|colistin
# 
# # Bray-Curtis dissimilarity test
# # Bray-Curtis dissimilarity test matrix
# dissimilarity_algorthmn <- "BrayCurtis_dissimilarity"
# bray_curtis_matrix <- vegdist(inPutData, method = "bray", na.rm = TRUE)
# 
# # Aitchison distance (also known as the Aitchison distance or CLR distance)
# # Aitchison distance test matrix
# dissimilarity_algorthmn <- "RobustAitchison" #"RobustAitchison"|"Aitchison"
# robust_aitchison_matrix <- vegdist(inPutData, method = "robust.aitchison", na.rm = TRUE)
# write.table(as.matrix(robust_aitchison_matrix), 
#             file = "./results_arg/r-analysis/data/robust_aitchison_matrix.tsv",
#             sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# 
# # Assign dissimilarity matrix
# dissimilarity_matrix <- bray_curtis_matrix #bray_curtis_matrix|robust_aitchison_matrix
# dissimilarity_algorthmn <- "BrayCurtis_dissimilarity" #"BrayCurtis_dissimilarity"|"RobustAitchison"|"Aitchison"
# 
# # Perform hierarchical clustering 
# hc <- hclust(dissimilarity_matrix, method = "ward.D2")
# dend <- as.dendrogram(hc, center = T)
# plot(dend, main = "Hierarchical Clustering Dendrogram", ylab = "Height",
#      xlab= paste(dissimilarity_algorthmn,"matrix",sep = " "))
# 
# # Ordination
# # Perform PCoA (cmdscale())
# dimentionality_reduction <- "PCoA"
# pcoa_results <- cmdscale(dissimilarity_matrix, eig = TRUE, k = 2)
# # Convert PCoA results to a dataframe
# pcoa_df0 <- as.data.frame(pcoa_results$points) %>% rename_all(~str_replace_all(., "V", "PCoA"))
# pcoa_df <- merge(rownames_to_column(pcoa_df0, var = "sampleID"),inMetadata)
# pcoA_Plot_base <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2)) + 
#   geom_point(aes(color = EstateOfOrigin, shape = socioeconomic_category), size = 3) + 
#   labs(title = "PCoA Plot", x = "PCoA Axis 1", y = "PCoA Axis 2") + 
#   theme(text = element_text(face = "plain", size = 20),
#         axis.text.y = element_text(angle = 0), legend.position="right")
# pcoA_Plot_base_ellipse <- pcoA_Plot_base +
#   stat_ellipse(aes(colour = socioeconomic_category), linewidth = 1) + # linewidth not size, since ggplot 3.4.0
#   coord_fixed(ratio = 1, clip = "off")
#   # linetype = socioeconomic_category,
# 
# # Calculate the species scores that represent the contribution of each species to the principal coordinates.
# # envfit function from the vegan package fits environmental variables (species) to the ordination
# pcoa_points <- pcoa_df0 %>% select(PCoA1, PCoA2) %>% as.matrix() %>% na.omit(.)
# # Remove columns with more than a certain percentage of NAs
# threshold <- 0.5  # Adjust this threshold as needed
# inPutData_filtered <- inPutData[, colMeans(is.na(inPutData)) < threshold]
# 
# # Impute remaining NAs
# inPutData_imputed <- setNames(
#   as.data.frame(lapply(inPutData_filtered, function(x) {
#     if (is.numeric(x)) {
#       x[is.na(x)] <- mean(x, na.rm = TRUE)
#     }
#     return(x)
#   })),
#   names(inPutData_filtered)
# )
# 
# envfit_results <- envfit(pcoa_results, inPutData_imputed, perm = 999, na.rm = TRUE)
# # Extract the scores from the envfit object and convert them to a data frame
# amrClass_scores <- as.data.frame(scores(envfit_results, display = "vectors")) %>%
#   rename_all(~str_replace_all(., "Dim", "PCoA"))
# amrClass_scores$amrClass <- rownames(amrClass_scores)
# # Calculate the length of vectors to determine importance 
# amrClass_scores$length <- sqrt(amrClass_scores$PCoA1^2 + amrClass_scores$PCoA2^2)
# amrClass_scores <- amrClass_scores %>% 
#   mutate(labelAngle = atan2(PCoA2, PCoA1) * 180 / pi, label_x = PCoA1 * 1.1, label_y = PCoA2 * 1.1)
# # Filter for the most important vectors (e.g., top 10) 
# top_n <- 100 # Adjust this number as needed 
# important_vectors <- amrClass_scores %>% arrange(desc(length)) %>% head(top_n)
# 
# # Calculate the scaling factor: PCoA1 and PCoA2 axes are much larger than the lengths of vectors
# # Causes vectors to appear very short and congested at the center.
# # Need to rescale the vectors to better match the scale of your PCoA axes
# division_fastor <- max(abs(c(pcoa_df$PCoA1, pcoa_df$PCoA2)))
# scaling_factor <- max(abs(c(important_vectors$PCoA1, important_vectors$PCoA2))) / division_fastor
# important_vectors$PCoA1_scaled <- important_vectors$PCoA1 / scaling_factor
# important_vectors$PCoA2_scaled <- important_vectors$PCoA2 / scaling_factor
# # Add the vectors to the plot: Use geom_segment to add the vectors and geom_text to label them
# pcoA_Plot <- pcoA_Plot_base_ellipse +
#   geom_segment(data = important_vectors, aes(x = 0, y = 0, xend = PCoA1_scaled, yend = PCoA2_scaled), 
#                arrow = arrow(length = unit(0.2, "cm")), color = "grey50") +
#   geom_label_repel(data = important_vectors, aes(x = PCoA1_scaled, y = PCoA2_scaled, label = amrClass, angle = labelAngle), 
#                    fill = "white", color = "black", box.padding = unit(0.3, "lines"), 
#                    point.padding = unit(0.3, "lines"), segment.color = "black", 
#                    size = 2) +
#   annotate("text", x = Inf, y = Inf, label = paste("Bi-plots rescaled by division by", scaling_factor),
#            hjust = 1, vjust = 1, size = 4, color = "grey50") + 
#   guides(color = guide_legend(override.aes = list(shape = 19, size = 3)), # Override shapes in legend
#          shape = guide_legend(override.aes = list(size = 3))) # Ensure shape sizes are consistent
# 
# # Display the plot
# pcoA_Plot
# 
# # Lets Save the plot now
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",dissimilarity_algorthmn,
#                "_important_vectors_Top",top_n,"_PCoA_Plot.png", sep = "")
# ggsave(fileName, plot = pcoA_Plot, width = 35, height = 25, units = "cm",limitsize = FALSE)
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",dissimilarity_algorthmn,
#                "_important_vectors_Top",top_n,"_PCoA_Plot.pdf", sep = "")
# ggsave(fileName, plot = pcoA_Plot, width = 35, height = 25, units = "cm",limitsize = FALSE,
#        device = cairo_pdf)
# 
# # Perform NMDS
# dimentionality_reduction <- "NMDS"
# nmds_results <- metaMDS(dissimilarity_matrix, k = 2, trymax = 100)
# # Convert PCoA results to a dataframe
# nmds_df0 <- as.data.frame(nmds_results$points) %>% rename_all(~str_replace_all(., "MDS", "NMDS"))
# nmds_df <- merge(rownames_to_column(nmds_df0, var = "sampleID"),inMetadata)
# nmds_Plot_base <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) + 
#   geom_point(aes(color = EstateOfOrigin, shape = socioeconomic_category), size = 3) + 
#   labs(title = "NMDS Plot", x = "NMDS Axis 1", y = "NMDS Axis 2") + 
#   theme(text = element_text(face = "plain", size = 20),
#         axis.text.y = element_text(angle = 0),
#         legend.position="right") 
# nmds_Plot_base_ellipse <- nmds_Plot_base +
#   stat_ellipse(aes(colour = socioeconomic_category), linewidth = 1) + # linewidth not size, since ggplot 3.4.0
#   coord_fixed(ratio = 1, clip = "off")
#   # linetype = socioeconomic_category, 
# # Calculate the species scores that represent the contribution of each species to the principal coordinates.
# # envfit function from the vegan package fits environmental variables (species) to the ordination
# nmds_points <- nmds_df0 %>% select(NMDS1, NMDS2) %>% as.matrix() %>% na.omit(.)
# 
# # Calculate environmental vectors
# envfit_results <- envfit(nmds_results, inPutData_imputed, perm = 999, na.rm = TRUE)
# # Extract the scores from the envfit object and convert them to a data frame
# amrClass_scores <- as.data.frame(scores(envfit_results, display = "vectors"))
# amrClass_scores$amrClass <- rownames(amrClass_scores)
# # Calculate the length of vectors to determine importance 
# amrClass_scores$length <- sqrt(amrClass_scores$NMDS1^2 + amrClass_scores$NMDS2^2)
# amrClass_scores <- amrClass_scores %>% 
#   mutate(angle = atan2(NMDS1, NMDS2) * 180 / pi, label_x = NMDS1 * 1.1, label_y = NMDS2 * 1.1)
# # Filter for the most important vectors (e.g., top 10) 
# top_n <- 100 # Adjust this number as needed 
# important_vectors <- amrClass_scores %>% arrange(desc(length)) %>% head(top_n)
# 
# # Calculate the scaling factor: PCoA1 and PCoA2 axes are much larger than the lengths of vectors
# # Causes vectors to appear very short and congested at the center.
# # Need to rescale the vectors to better match the scale of your PCoA axes
# division_fastor <- max(abs(c(nmds_df$NMDS1, nmds_df$NMDS2)))
# scaling_factor <- max(abs(c(important_vectors$NMDS1, important_vectors$NMDS2))) / division_fastor
# important_vectors$NMDS1_scaled <- important_vectors$NMDS1 / scaling_factor
# important_vectors$NMDS2_scaled <- important_vectors$NMDS2 / scaling_factor
# 
# # Add the vectors to the plot: Use geom_segment to add the vectors and geom_text to label them
# nmds_Plot <- nmds_Plot_base_ellipse +
#   geom_segment(data = important_vectors, aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled), 
#                arrow = arrow(length = unit(0.2, "cm")), color = "grey50") +
#   geom_label_repel(data = important_vectors, aes(x = NMDS1_scaled, y = NMDS2_scaled, label = amrClass), 
#                    fill = "white", color = "black", box.padding = unit(0.3, "lines"), 
#                    point.padding = unit(0.3, "lines"), segment.color = "black", 
#                    size = 2, max.overlaps = 30) +
#   annotate("text", x = Inf, y = Inf, label = paste("Bi-plots rescaled by division by", scaling_factor),
#            hjust = 1, vjust = 1, size = 4, color = "grey50")
# 
# # Display the plot
# nmds_Plot
# 
# # Lets Save the plot now
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",
#                dissimilarity_algorthmn,"_important_vectors_Top",top_n,"_NMDS_Plot.png", sep = "")
# ggsave(fileName, plot = nmds_Plot, width = 35, height = 25, units = "cm",limitsize = FALSE)
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",
#                dissimilarity_algorthmn,"_important_vectors_Top",top_n,"_NMDS_Plot.pdf", sep = "")
# ggsave(fileName, plot = nmds_Plot, width = 35, height = 25, units = "cm",limitsize = FALSE)
# 
# 
# # PERMANOVA
# # Performing PERMANOVA (Permutational Multivariate Analysis of Variance) is a powerful way to test for 
# # significant differences in community composition among groups. Helps explain if different groups 
# # (e.g., habitats, treatments, or time points) have significantly different compositions of species or other multivariate data.
# # Assuming 'socioeconomic_category' is a factor indicating the group (low_income, middle_income, high_income)
# permanova_results <- adonis2(dissimilarity_matrix ~ socioeconomic_category + EstateOfOrigin + weekNo + Seqrun, 
#                              data = inMetadata, permutations = 999, by = "margin")
# print(permanova_results)
# 
# # Save the results to a text file
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",
#                dissimilarity_algorthmn,"_PERMANOVA_results.txt", sep = "")
# writeLines(capture.output((permanova_results)), con = fileName, sep = "\n", useBytes = TRUE)
# 
# # Assuming 'socioeconomic_category' is a factor indicating the group (low_income, middle_income, high_income)
# permanova_results0 <- adonis2(dissimilarity_matrix ~ socioeconomic_category + weekNo + Seqrun, 
#                              data = inMetadata, permutations = 999, by = "margin")
# print(permanova_results0)
# 
# # Save the results to a text file
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",
#                dissimilarity_algorthmn,"_PERMANOVA_margin.txt", sep = "")
# writeLines(capture.output((permanova_results0)), con = fileName, sep = "\n", useBytes = TRUE)
# 
# # Assuming 'socioeconomic_category' is a factor indicating the group (low_income, middle_income, high_income)
# permanova_results1 <- adonis2(dissimilarity_matrix ~ socioeconomic_category + EstateOfOrigin + weekNo + Seqrun, 
#                              data = inMetadata, permutations = 999, by = "terms")
# print(permanova_results1)
# 
# # Save the results to a text file
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",
#                dissimilarity_algorthmn,"_PERMANOVATerms_results.txt", sep = "")
# writeLines(capture.output((permanova_results1)), con = fileName, sep = "\n", useBytes = TRUE)
# 
# 
# # Assuming 'socioeconomic_category' is a factor indicating the group (low_income, middle_income, high_income)
# permanova_results2 <- adonis2(dissimilarity_matrix ~ socioeconomic_category, 
#                              data = inMetadata, permutations = 999)
# print(permanova_results2)
# 
# # Save the results to a text file
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",
#                dissimilarity_algorthmn,"_PERMANOVA_SCWk_results.txt", sep = "")
# writeLines(capture.output((permanova_results2)), con = fileName, sep = "\n", useBytes = TRUE)
# 
# # Assuming 'socioeconomic_category' is a factor indicating the group (low_income, middle_income, high_income)
# # Using week as a strata, to test the effect of distance and socio-economic class on ---
# # your response variable while accounting for the variability within each week.
# permanova_results3 <- adonis2(dissimilarity_matrix ~ socioeconomic_category, 
#                               data = inMetadata, permutations = 999,
#                               strata = inMetadata$weekNo)
# print(permanova_results3)
# 
# # Save the results to a text file
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",
#                dissimilarity_algorthmn,"_PERMANOVA_strataweekNo_results.txt", sep = "")
# writeLines(capture.output((permanova_results3)), con = fileName, sep = "\n", useBytes = TRUE)
# 
# # Pairwise adonis : pairwiseAdonis package
# # Assuming 'socioeconomic_category' is a factor indicating the group (low_income, middle_income, high_income)
# pairwiseAdonis_results <- pairwise.adonis(dissimilarity_matrix, inMetadata$socioeconomic_category)
# print(pairwiseAdonis_results)
# 
# # Save the results to a text file
# fileName=paste("./results_arg/r-analysis/plots/",shortDate,"_",inPutDateName,"_",
#                dissimilarity_algorthmn,"_pairwiseAdonis_results.txt", sep = "")
# writeLines(capture.output((pairwiseAdonis_results)), con = fileName, sep = "\n", useBytes = TRUE)
