## =============================================================================
## 06_ARG_Ordination_PERMANOVA_AllMetrics.R
## -----------------------------------------------------------------------------
## Runs the full ordination + PERMANOVA workflow (Bray-Curtis & Robust
## Aitchison distances -> hierarchical clustering -> PCoA -> NMDS ->
## PERMANOVA suite -> pairwise PERMANOVA) across EVERY combination of:
##   - abundance metric:  relativeAbundance | RPKM | RPK
##   - taxonomic level:   DrugClass | AMRClass | AMRGeneFamily
##
## This replaces the earlier single-matrix / per-drug-class scripts (which
## duplicated ~300 lines of PCoA/NMDS/PERMANOVA code per dataset) with one
## parametrised loop. Adding a new metric or level requires no code changes --
## only editing the `abundance_metrics` / `taxonomic_levels` vectors below.
##
## See companion doc: ANALYSIS_METHODS_AND_INTERPRETATION.md for how to read
## and report every output this script produces.
## =============================================================================

rm(list = ls())
unlink(".Rhistory")
gc()
setwd("/home/gkibet/hpc_mnt/bioinformatics/github/metagenomics/data/visualization/")
lib <- "/home/gkibet/R/x86_64-pc-linux-gnu-library/4.4"
source("/home/gkibet/hpc_mnt/bioinformatics/github/metagenomics/data/visualization/scripts/functions.R")
source("/home/gkibet/hpc_mnt/bioinformatics/github/metagenomics/data/visualization/scripts/00.ordination_permanova_functions.R")

requiredCRANPackages <- c(
  "dplyr", "tidyr", "tibble", "ggplot2", "stringr", "vegan", "ggrepel",
  "dendextend", "pairwiseAdonis", "readr")
requiredGITHUBPackages <- c("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
installNloadCRANpackages(requiredPackages = requiredCRANPackages, lib = lib)
# Uncomment on first run only, or if pairwiseAdonis isn't yet installed:
# installNloadGitHubpackages(requiredPackages = requiredGITHUBPackages, lib = lib)

shortDate <- "20260716"

# ---------------------------------------------------------------------------
# Parameters -- edit these to control what gets analysed
# ---------------------------------------------------------------------------
resfinderdatabasePath <- "./results_arg/r-analysis/data/"
plotsbasePath          <- "./results_arg/r-analysis/plots/"

# NOTE on "RPK": true RPK (reads-per-kilobase, no library-size / depth
# correction) is NOT normalised for sequencing depth across samples.
# It is included here for completeness per your request, but cross-sample
# comparisons (ordination, PERMANOVA) on RPK_* matrices should be interpreted
# cautiously -- RPKM and relative abundance are the depth-normalised metrics
# intended for this kind of between-sample comparison. See the companion
# methods document for the full rationale.
abundance_metrics <- c("relativeAbundance", "RPKM", "RPK")
taxonomic_levels  <- c("DrugClass", "AMRClass", "AMRGeneFamily")

# Non-numeric metadata columns that may be embedded in the abundance matrix
# files and must be dropped before distance/ordination calculations
nonNumColNames <- c("EstateOfOrigin", "year_Week", "CountyOfOrigin",
                    "weekNo", "weekNoPerCountyOfOrigin", "Seqrun")

na_threshold  <- 0.5   # drop amrClass/drugClass columns with >50% missing before envfit
top_n_vectors <- 100   # number of envfit vectors labelled on each biplot

distance_methods <- c(BrayCurtis_dissimilarity = "bray",
                      RobustAitchison          = "robust.aitchison")

# ---------------------------------------------------------------------------
# Metadata (shared across all metric x level combinations)
# ---------------------------------------------------------------------------
inMetadata_file <- paste0(resfinderdatabasePath, shortDate, "_AMRData_Nmetadata.tsv")
cat("Reading metadata:", inMetadata_file, "\n")
inMetadata <- read.csv(inMetadata_file, sep = "\t", header = TRUE)

# ---------------------------------------------------------------------------
# Master log -- all console/model output across every combination is appended
# here, in addition to the per-model .txt files saved alongside each plot
# ---------------------------------------------------------------------------
logFileName <- paste0("./results_arg/r-analysis/OrdinationPERMANOVA_AllMetrics_",
                      format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
writeLines("", logFileName)

# Consolidated NMDS stress summary across every metric x level x distance
# combination -- saved once at the end so you don't have to open every plot
# to check whether each NMDS representation was trustworthy (stress < 0.2)
stress_summary <- data.frame(
  abundance_metric = character(), taxonomic_level = character(),
  distance_method = character(), nmds_stress = numeric(),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# Main loop: abundance metric x taxonomic level
# ---------------------------------------------------------------------------
for (metric in abundance_metrics) {
  for (level in taxonomic_levels) {
    
    inPutDateName <- paste0(metric, "Per", level, "_mtx")
    matrix_file <- paste0(resfinderdatabasePath, shortDate, "_AMRData_", metric, "Per", level, "_Nmetadata.tsv")
    
    if (!file.exists(matrix_file)) {
      cat("SKIP:", matrix_file, "not found -- skipping", metric, "x", level, "\n")
      next
    }
    
    cat("=======================================================\n")
    cat("Processing:", metric, "x", level, "\n")
    cat("=======================================================\n")
    
    sink(file = logFileName, append = TRUE)
    
    out_dir <- paste0(plotsbasePath, metric, "_", level, "/")
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    
    # ---- Read + clean the abundance matrix ----
    inPutData00 <- read.csv(matrix_file, sep = "\t", header = TRUE) %>%
      tibble::column_to_rownames(var = "sampleID")
    inPutData <- inPutData00 %>% dplyr::select(-dplyr::any_of(nonNumColNames))
    
    inMetadata_sub <- inMetadata %>% dplyr::filter(sampleID %in% rownames(inPutData))
    
    # ---- NA filtering + imputation (shared input for both PCoA & NMDS envfit) ----
    prepped <- filter_and_impute(inPutData, na_threshold = na_threshold)
    inPutData_filtered <- prepped$filtered
    inPutData_imputed  <- prepped$imputed
    
    # ---- Distance matrices ----
    # ---- Distance matrices ----
    inPutData_numeric <- inPutData %>% dplyr::select(-(1:2))
    dist_matrices <- build_distance_matrices(inPutData_numeric, methods = distance_methods)
    for (nm in names(dist_matrices)) {
      write.table(as.matrix(dist_matrices[[nm]]),
                  file = paste0(out_dir, shortDate, "_", nm, "_", inPutDateName, ".tsv"),
                  sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
    }
    
    # ---- Per-distance-matrix analysis: dendrogram, PCoA, NMDS, PERMANOVA ----
    for (dissimilarity_algorthmn in names(dist_matrices)) {
      cat("  Distance method:", dissimilarity_algorthmn, "\n")
      dissimilarity_matrix <- dist_matrices[[dissimilarity_algorthmn]]
      out_prefix <- paste0(out_dir, shortDate, "_", inPutDateName, "_", dissimilarity_algorthmn)
      
      # Hierarchical clustering dendrogram
      hc <- hclust(dissimilarity_matrix, method = "ward.D2")
      png(paste0(out_prefix, "_dendrogram.png"), width = 25, height = 15, units = "cm", res = 300)
      plot(as.dendrogram(hc, center = TRUE), main = "Hierarchical Clustering Dendrogram",
           ylab = "Height", xlab = paste(dissimilarity_algorthmn, "matrix"))
      dev.off()
      
      # PCoA + envfit biplot
      cat("    PCoA ordination + envfit biplot\n")
      pcoa_out <- build_ordination_biplot("PCoA", dissimilarity_matrix, inMetadata_sub,
                                          inPutData_imputed, top_n = top_n_vectors)
      ggsave(paste0(out_prefix, "_PCoA_Plot.png"), plot = pcoa_out$plot,
             width = 35, height = 25, units = "cm", limitsize = FALSE)
      ggsave(paste0(out_prefix, "_PCoA_Plot.pdf"), plot = pcoa_out$plot,
             width = 35, height = 25, units = "cm", limitsize = FALSE, device = cairo_pdf)
      
      # NMDS + envfit biplot
      cat("    NMDS ordination + envfit biplot\n")
      nmds_out <- build_ordination_biplot("NMDS", dissimilarity_matrix, inMetadata_sub,
                                          inPutData_imputed, top_n = top_n_vectors)
      ggsave(paste0(out_prefix, "_NMDS_Plot.png"), plot = nmds_out$plot,
             width = 35, height = 25, units = "cm", limitsize = FALSE)
      ggsave(paste0(out_prefix, "_NMDS_Plot.pdf"), plot = nmds_out$plot,
             width = 35, height = 25, units = "cm", limitsize = FALSE)
      
      stress_summary <- rbind(stress_summary, data.frame(
        abundance_metric = metric, taxonomic_level = level,
        distance_method = dissimilarity_algorthmn, nmds_stress = nmds_out$stress
      ))
      
      # PERMANOVA suite (5 models) + pairwise adonis
      cat("    PERMANOVA suite\n")
      permanova_out <- run_permanova_suite(dissimilarity_matrix, inMetadata_sub, out_prefix)
    }
    
    sink()
  }
}

cat("All abundance-metric x taxonomic-level ordination/PERMANOVA analyses complete.\n")
cat("See:", logFileName, "for the full combined run log.\n")

# Save the consolidated NMDS stress summary -- check this before trusting
# any NMDS plot; rows with nmds_stress >= 0.2 should be flagged/caveated
stress_summary_file <- paste0(plotsbasePath, shortDate, "_NMDS_stress_summary.tsv")
write.table(stress_summary, file = stress_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("NMDS stress summary saved to:", stress_summary_file, "\n")
print(stress_summary)
