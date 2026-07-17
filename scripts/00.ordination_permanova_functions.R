## =============================================================================
## ordination_permanova_functions.R
## -----------------------------------------------------------------------------
## Reusable helper functions for AMR abundance ordination + PERMANOVA testing.
## Source this file from any driver script with:
##   source(".../ordination_permanova_functions.R")
##
## These functions replace the copy-pasted PCoA/NMDS/PERMANOVA blocks that were
## duplicated across the original per-matrix and per-drug-class scripts, so a
## single change here (e.g. plot styling, a new PERMANOVA model) propagates
## everywhere it's used.
## =============================================================================

# ---------------------------------------------------------------------------
# 1. Build one or more vegan::vegdist() dissimilarity/distance matrices
# ---------------------------------------------------------------------------
# methods: a NAMED vector, names = label used in filenames/output,
#          values = the vegdist() `method` argument
build_distance_matrices <- function(mtx,
                                    methods = c(BrayCurtis_dissimilarity = "bray",
                                                RobustAitchison          = "robust.aitchison")) {
  out <- list()
  for (nm in names(methods)) {
    cat("    Building", nm, "distance matrix (method =", methods[[nm]], ")...\n")
    out[[nm]] <- vegan::vegdist(mtx, method = methods[[nm]], na.rm = TRUE)
  }
  out
}

# ---------------------------------------------------------------------------
# 2. Shared NA-column filtering + mean imputation
#    (used to build the envfit() input matrix for both PCoA and NMDS biplots)
# ---------------------------------------------------------------------------
filter_and_impute <- function(mtx, na_threshold = 0.5) {
  mtx_filtered <- mtx[, colMeans(is.na(mtx)) < na_threshold, drop = FALSE]
  mtx_imputed <- setNames(
    as.data.frame(lapply(mtx_filtered, function(x) {
      if (is.numeric(x)) x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    })),
    names(mtx_filtered)
  )
  list(filtered = mtx_filtered, imputed = mtx_imputed)
}

# ---------------------------------------------------------------------------
# 3. Generic ordination + envfit biplot builder
#    ordination_type = "PCoA" (cmdscale) or "NMDS" (metaMDS)
#    Produces the same style of plot used in the original scripts: point plot
#    colored/shaped by metadata, stat_ellipse by group, envfit vectors for the
#    top_n most important AMR classes, rescaled to fit the ordination axes.
# ---------------------------------------------------------------------------
build_ordination_biplot <- function(ordination_type = c("PCoA", "NMDS"),
                                    dissimilarity_matrix,
                                    inMetadata,
                                    inPutData_imputed,
                                    top_n = 100,
                                    colour_var = "EstateOfOrigin",
                                    shape_var = "socioeconomic_category",
                                    ellipse_var = "socioeconomic_category",
                                    feature_label = "amrClass") {
  
  ordination_type <- match.arg(ordination_type)
  
  if (ordination_type == "PCoA") {
    ord_results <- cmdscale(dissimilarity_matrix, eig = TRUE, k = 2)
    axis_prefix <- "PCoA"
    df0 <- as.data.frame(ord_results$points) %>%
      dplyr::rename_all(~stringr::str_replace_all(., "V", axis_prefix))
    envfit_dim_prefix <- "Dim"
    stress_value <- NA_real_   # not applicable to PCoA
  } else {
    ord_results <- vegan::metaMDS(dissimilarity_matrix, k = 2, trymax = 100)
    axis_prefix <- "NMDS"
    df0 <- as.data.frame(ord_results$points) %>%
      dplyr::rename_all(~stringr::str_replace_all(., "MDS", axis_prefix))
    envfit_dim_prefix <- "NMDS"
    stress_value <- ord_results$stress
    cat("      NMDS stress =", round(stress_value, 4),
        ifelse(stress_value < 0.2, "(acceptable, <0.2)", "(CAUTION: >=0.2, interpret with care)"), "\n")
  }
  
  ax1 <- paste0(axis_prefix, "1")
  ax2 <- paste0(axis_prefix, "2")
  df <- merge(tibble::rownames_to_column(df0, var = "sampleID"), inMetadata)
  
  base_plot <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[ax1]], y = .data[[ax2]])) +
    ggplot2::geom_point(ggplot2::aes(color = .data[[colour_var]], shape = .data[[shape_var]]), size = 3) +
    ggplot2::labs(title = paste(ordination_type, "Plot"),
                  x = paste(axis_prefix, "Axis 1"), y = paste(axis_prefix, "Axis 2")) +
    ggplot2::theme(text = ggplot2::element_text(face = "plain", size = 20),
                   axis.text.y = ggplot2::element_text(angle = 0), legend.position = "right")
  
  base_plot_ellipse <- base_plot +
    ggplot2::stat_ellipse(ggplot2::aes(colour = .data[[ellipse_var]]), linewidth = 1) +
    ggplot2::coord_fixed(ratio = 1, clip = "off")
  
  # envfit: fits each AMR class as a "species vector" onto the ordination axes,
  # i.e. which classes' abundance gradients align most with each axis
  envfit_results <- vegan::envfit(ord_results, inPutData_imputed, perm = 999, na.rm = TRUE)
  scores_df <- as.data.frame(vegan::scores(envfit_results, display = "vectors")) %>%
    dplyr::rename_all(~stringr::str_replace_all(., envfit_dim_prefix, axis_prefix))
  scores_df[[feature_label]] <- rownames(scores_df)
  scores_df$length <- sqrt(scores_df[[ax1]]^2 + scores_df[[ax2]]^2)
  scores_df <- scores_df %>%
    dplyr::mutate(labelAngle = atan2(.data[[ax2]], .data[[ax1]]) * 180 / pi)
  
  important_vectors <- scores_df %>% dplyr::arrange(dplyr::desc(length)) %>% head(top_n)
  
  # Rescale vectors so they're visible relative to the (much larger) sample-point spread
  division_factor <- max(abs(c(df[[ax1]], df[[ax2]])))
  scaling_factor <- max(abs(c(important_vectors[[ax1]], important_vectors[[ax2]]))) / division_factor
  important_vectors[[paste0(ax1, "_scaled")]] <- important_vectors[[ax1]] / scaling_factor
  important_vectors[[paste0(ax2, "_scaled")]] <- important_vectors[[ax2]] / scaling_factor
  
  final_plot <- base_plot_ellipse +
    ggplot2::geom_segment(
      data = important_vectors,
      ggplot2::aes(x = 0, y = 0,
                   xend = .data[[paste0(ax1, "_scaled")]],
                   yend = .data[[paste0(ax2, "_scaled")]]),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")), color = "grey50") +
    ggrepel::geom_label_repel(
      data = important_vectors,
      ggplot2::aes(x = .data[[paste0(ax1, "_scaled")]], y = .data[[paste0(ax2, "_scaled")]],
                   label = .data[[feature_label]]),
      fill = "white", color = "black", box.padding = ggplot2::unit(0.3, "lines"),
      point.padding = ggplot2::unit(0.3, "lines"), segment.color = "black",
      size = 2, max.overlaps = 30) +
    ggplot2::annotate("text", x = Inf, y = Inf,
                      label = paste("Bi-plot vectors rescaled by division by", round(scaling_factor, 3)),
                      hjust = 1, vjust = 1, size = 4, color = "grey50") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 19, size = 3)),
                    shape = ggplot2::guide_legend(override.aes = list(size = 3)))
  
  # For NMDS, print the stress value directly on the plot -- this is the
  # number a reader needs to judge whether the 2D representation is trustworthy
  if (ordination_type == "NMDS") {
    final_plot <- final_plot +
      ggplot2::annotate("text", x = Inf, y = -Inf,
                        label = paste0("Stress = ", round(stress_value, 4),
                                       ifelse(stress_value < 0.2, "", "  (CAUTION: >=0.2)")),
                        hjust = 1, vjust = -0.5, size = 4,
                        color = ifelse(stress_value < 0.2, "grey30", "firebrick"))
  }
  
  list(ordination_object = ord_results, coordinates_df = df, envfit = envfit_results,
       important_vectors = important_vectors, plot = final_plot, scaling_factor = scaling_factor,
       stress = stress_value)
}

# ---------------------------------------------------------------------------
# 4. PERMANOVA test suite -- runs the same battery of models used previously:
#    - full model, marginal (type III-like) test
#    - model excluding EstateOfOrigin (avoids confound w/ socioeconomic_category)
#    - full model, sequential (type I) test
#    - socioeconomic_category alone
#    - socioeconomic_category with weekNo as permutation strata
#    - pairwise adonis between socioeconomic_category levels
# ---------------------------------------------------------------------------
run_permanova_suite <- function(dissimilarity_matrix, inMetadata, out_prefix, n_perm = 999) {
  results <- list()
  
  cat("    PERMANOVA: full model ~ margin\n")
  results$full_margin <- vegan::adonis2(
    dissimilarity_matrix ~ socioeconomic_category + EstateOfOrigin + weekNo + Seqrun,
    data = inMetadata, permutations = n_perm, by = "margin")
  writeLines(capture.output(results$full_margin), paste0(out_prefix, "_PERMANOVA_full_margin.txt"))
  
  cat("    PERMANOVA: no-estate model ~ margin\n")
  results$noEstate_margin <- vegan::adonis2(
    dissimilarity_matrix ~ socioeconomic_category + weekNo + Seqrun,
    data = inMetadata, permutations = n_perm, by = "margin")
  writeLines(capture.output(results$noEstate_margin), paste0(out_prefix, "_PERMANOVA_noEstate_margin.txt"))
  
  cat("    PERMANOVA: full model ~ terms (sequential)\n")
  results$full_terms <- vegan::adonis2(
    dissimilarity_matrix ~ socioeconomic_category + EstateOfOrigin + weekNo + Seqrun,
    data = inMetadata, permutations = n_perm, by = "terms")
  writeLines(capture.output(results$full_terms), paste0(out_prefix, "_PERMANOVA_terms.txt"))
  
  cat("    PERMANOVA: socioeconomic_category only\n")
  results$socio_only <- vegan::adonis2(
    dissimilarity_matrix ~ socioeconomic_category,
    data = inMetadata, permutations = n_perm)
  writeLines(capture.output(results$socio_only), paste0(out_prefix, "_PERMANOVA_socioOnly.txt"))
  
  cat("    PERMANOVA: socioeconomic_category, strata = weekNo\n")
  results$socio_strataWeek <- vegan::adonis2(
    dissimilarity_matrix ~ socioeconomic_category,
    data = inMetadata, permutations = n_perm, strata = inMetadata$weekNo)
  writeLines(capture.output(results$socio_strataWeek), paste0(out_prefix, "_PERMANOVA_strataWeek.txt"))
  
  cat("    Pairwise adonis: socioeconomic_category\n")
  results$pairwise <- pairwiseAdonis::pairwise.adonis(dissimilarity_matrix, inMetadata$socioeconomic_category)
  writeLines(capture.output(results$pairwise), paste0(out_prefix, "_pairwiseAdonis.txt"))
  
  results
}