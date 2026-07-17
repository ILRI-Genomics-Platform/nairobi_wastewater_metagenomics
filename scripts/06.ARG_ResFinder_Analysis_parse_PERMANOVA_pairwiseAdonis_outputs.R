## =============================================================================
## 09_parse_PERMANOVA_pairwiseAdonis_outputs.R
## -----------------------------------------------------------------------------
## Consolidates every PERMANOVA (adonis2) and pairwiseAdonis .txt output file
## produced by script 07 into TWO tidy summary tables:
##   1. permanova_summary.tsv   -- one row per (metric, level, distance, model, term)
##   2. pairwiseAdonis_summary.tsv -- one row per (metric, level, distance, pair)
##
## Filenames are parsed directly, so this expects the naming convention used
## by 07_ARG_Ordination_PERMANOVA_AllMetrics.R:
##   {shortDate}_{metric}Per{level}_mtx_{distance}_PERMANOVA_{model}.txt
##   {shortDate}_{metric}Per{level}_mtx_{distance}_pairwiseAdonis.txt
## =============================================================================

library(dplyr)
library(stringr)
library(purrr)
library(readr)

# ---------------------------------------------------------------------------
# Parameters -- point this at wherever your .txt outputs actually live
# ---------------------------------------------------------------------------
results_dir <- "./results_arg/r-analysis/plots/"   # searched recursively
out_dir     <- "./results_arg/r-analysis/data/"

# ---------------------------------------------------------------------------
# Helper: parse metric / level / distance out of a filename
# ---------------------------------------------------------------------------
parse_filename <- function(filename) {
  base <- basename(filename)
  # Strip the leading date and trailing "_mtx_..." to isolate metric+level
  m <- str_match(base, "^\\d+_([A-Za-z]+)Per([A-Za-z]+)_mtx_([A-Za-z_]+?)_(PERMANOVA_.*|pairwiseAdonis)\\.txt$")
  if (is.na(m[1, 1])) return(NULL)
  list(metric = m[1, 2], level = m[1, 3], distance = m[1, 4], suffix = m[1, 5])
}

# ---------------------------------------------------------------------------
# Parse ONE adonis2 PERMANOVA .txt file into a tidy data frame (one row/term)
# ---------------------------------------------------------------------------
parse_permanova_file <- function(filepath) {
  meta <- parse_filename(filepath)
  if (is.null(meta)) {
    cat("  Skipping (name doesn't match expected pattern):", basename(filepath), "\n")
    return(NULL)
  }
  model_name <- str_remove(meta$suffix, "^PERMANOVA_")
  
  lines <- readLines(filepath, warn = FALSE)
  # The results table sits between the "adonis2(formula" line and the "---" line
  start <- which(str_detect(lines, "^vegan::adonis2")) + 1
  end_candidates <- which(str_detect(lines, "^---"))
  end <- if (length(end_candidates) > 0) min(end_candidates) - 1 else length(lines)
  if (length(start) == 0 || start > end) return(NULL)
  
  table_lines <- lines[start:end]
  table_lines <- table_lines[nchar(trimws(table_lines)) > 0]
  
  # First remaining line is the header; parse fixed-width-ish whitespace split
  header <- str_split(str_trim(table_lines[1]), "\\s+")[[1]]
  data_lines <- table_lines[-1]
  
  rows <- map(data_lines, function(l) {
    # Term name may itself contain no spaces (matches your factor names), so
    # splitting on whitespace after trimming works reliably here
    parts <- str_split(str_trim(l), "\\s+")[[1]]
    term <- parts[1]
    values <- suppressWarnings(as.numeric(parts[2:length(parts)]))
    # adonis2 columns are: Df, SumOfSqs, R2, F, Pr(>F) -- but the last two
    # (F, Pr(>F)) are absent for Residual/Total rows
    tibble(
      term = term,
      Df = values[1],
      SumOfSqs = values[2],
      R2 = values[3],
      F_model = if (length(values) >= 4) values[4] else NA_real_,
      p_value = if (length(values) >= 5) values[5] else NA_real_
    )
  })
  
  bind_rows(rows) %>%
    mutate(abundance_metric = meta$metric, taxonomic_level = meta$level,
           distance_method = meta$distance, permanova_model = model_name,
           source_file = basename(filepath)) %>%
    relocate(abundance_metric, taxonomic_level, distance_method, permanova_model)
}

# ---------------------------------------------------------------------------
# Parse ONE pairwiseAdonis .txt file into a tidy data frame (one row/pair)
#
# IMPORTANT: R's capture.output()/print() wraps wide data frames into
# MULTIPLE blocks when they exceed console width. In practice this means
# `p.adjusted` and `sig` live in a SEPARATE block further down the file,
# not on the same line as pair/Df/SumsOfSqs/F.Model/R2/p.value. E.g.:
#
#                          pairs Df SumsOfSqs   F.Model         R2 p.value
#  1         high_income vs mixed  1 0.2189116  1.940959 0.04984364   0.121
#  ...
#    p.adjusted sig
#  1      0.726
#  2      0.006   *
#
# This function detects block boundaries (a header line = first token is
# NOT parseable as an integer row index) and re-joins wrapped blocks back
# together by row index.
# ---------------------------------------------------------------------------
parse_pairwise_file <- function(filepath) {
  meta <- parse_filename(filepath)
  if (is.null(meta)) {
    cat("  Skipping (name doesn't match expected pattern):", basename(filepath), "\n")
    return(NULL)
  }
  
  lines <- readLines(filepath, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]
  
  is_header_line <- function(l) {
    first_tok <- str_split(str_trim(l), "\\s+")[[1]][1]
    is.na(suppressWarnings(as.integer(first_tok)))
  }
  header_idx <- which(vapply(lines, is_header_line, logical(1)))
  if (length(header_idx) == 0) return(NULL)
  block_bounds <- c(header_idx, length(lines) + 1)
  
  # ---- First block: row_id, pair (contains " vs "), Df, SumsOfSqs, F.Model, R2, p.value ----
  first_data <- lines[(header_idx[1] + 1):(block_bounds[2] - 1)]
  first_rows <- map(first_data, function(l) {
    parts <- str_split(str_trim(l), "\\s+")[[1]]
    row_id <- as.integer(parts[1])
    vs_idx <- which(parts == "vs")
    if (length(vs_idx) == 0) return(NULL)
    pair_label <- paste(parts[2:(vs_idx + 1)], collapse = " ")
    remainder <- suppressWarnings(as.numeric(parts[(vs_idx + 2):length(parts)]))
    tibble(row_id = row_id, pair = pair_label,
           Df = remainder[1], SumsOfSqs = remainder[2],
           F_model = remainder[3], R2 = remainder[4], p_value = remainder[5])
  })
  result <- bind_rows(first_rows)
  
  # ---- Any additional wrapped blocks (e.g. p.adjusted, sig) -- join by row_id ----
  if (length(header_idx) > 1) {
    for (b in 2:length(header_idx)) {
      header_tokens <- str_split(str_trim(lines[header_idx[b]]), "\\s+")[[1]]
      data_lines <- lines[(header_idx[b] + 1):(block_bounds[b + 1] - 1)]
      
      blk_rows <- map(data_lines, function(l) {
        parts <- str_split(str_trim(l), "\\s+")[[1]]
        row_id <- as.integer(parts[1])
        values <- parts[-1]   # remaining tokens, may be fewer than header_tokens
        # (e.g. "sig" column blank for non-significant rows)
        row <- setNames(as.list(values), header_tokens[seq_along(values)])
        as_tibble(c(list(row_id = row_id), row))
      })
      blk_df <- bind_rows(blk_rows)
      # numeric columns (p.adjusted etc.) parse to numeric; leave "sig" (asterisks) as character
      num_cols <- setdiff(names(blk_df), c("row_id", "sig"))
      blk_df <- blk_df %>% mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.))))
      
      result <- dplyr::left_join(result, blk_df, by = "row_id")
    }
  }
  
  result %>%
    dplyr::select(-row_id) %>%
    dplyr::rename_with(~ "p_adjusted", any_of("p.adjusted")) %>%
    mutate(abundance_metric = meta$metric, taxonomic_level = meta$level,
           distance_method = meta$distance, source_file = basename(filepath)) %>%
    relocate(abundance_metric, taxonomic_level, distance_method)
}

# ---------------------------------------------------------------------------
# Run across every matching file in results_dir (recursive)
# ---------------------------------------------------------------------------
permanova_files <- list.files(results_dir, pattern = "_PERMANOVA_.*\\.txt$",
                              recursive = TRUE, full.names = TRUE)
pairwise_files  <- list.files(results_dir, pattern = "_pairwiseAdonis\\.txt$",
                              recursive = TRUE, full.names = TRUE)

cat("Found", length(permanova_files), "PERMANOVA files and",
    length(pairwise_files), "pairwiseAdonis files\n")

permanova_summary <- map_dfr(permanova_files, parse_permanova_file)
pairwise_summary  <- map_dfr(pairwise_files, parse_pairwise_file)

# ---------------------------------------------------------------------------
# Save consolidated tables
# ---------------------------------------------------------------------------
write_tsv(permanova_summary, paste0(out_dir, "PERMANOVA_summary_ALL.tsv"))
write_tsv(pairwise_summary, paste0(out_dir, "pairwiseAdonis_summary_ALL.tsv"))

cat("Saved:\n  ", paste0(out_dir, "PERMANOVA_summary_ALL.tsv"), "\n  ",
    paste0(out_dir, "pairwiseAdonis_summary_ALL.tsv"), "\n")

# ---------------------------------------------------------------------------
# Quick sanity views -- helpful right after running this
# ---------------------------------------------------------------------------
cat("\n--- socioeconomic_category rows only, noEstate_margin & socioOnly models ---\n")
permanova_summary %>%
  filter(term == "socioeconomic_category",
         permanova_model %in% c("noEstate_margin", "socioOnly")) %>%
  arrange(taxonomic_level, abundance_metric, distance_method) %>%
  print(n = Inf)

cat("\n--- Pairs significant at p.adjusted < 0.05, any metric/level/distance ---\n")
pairwise_summary %>%
  filter(p_adjusted < 0.05) %>%
  arrange(taxonomic_level, abundance_metric, distance_method, p_adjusted) %>%
  print(n = Inf)
