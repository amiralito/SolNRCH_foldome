#!/usr/bin/env Rscript
################################################################################
# NRC Heatmap Analysis
#
# Author: AmirAli Toghani
# Repository: https://github.com/amiralito/SolNRCH_foldome
#
# Generates the per-clade SNI heatmaps published in the paper:
#   - Transposed pheatmap (metrics on rows, clades on columns)
#   - annotation_row for metric source (Scientist / Both / Co-Scientist) and
#     structural category; annotation_col for clade n_members
#   - Blue-white-red palette (#1C39BB / white / #ED2939), breaks -3 to 3
#   - cellwidth=30, cellheight=30, fontsize=11, angle_col=45, border_color="white"
#   - Two heatmaps per dataset: full set (all metrics) and parameter subsets
#     (Scientist + overlap, Co-Scientist + overlap)
#
# Inputs (relative to this script):
#   main/resistosome_analysis_summary.xlsx
#   test/resistosome_analysis_summary.xlsx
#   ../phylo/clades/*.tree
#
# Outputs (PDF, PNG, SVG):
#   main/heatmaps/    main_clade_heatmap{,_zraw}.{pdf,png,svg}
#   main/comparison/  main_{scientist,coscientist}_params.{pdf,png,svg}
#   test/heatmaps/    test_clade_heatmap{,_zraw}.{pdf,png,svg}
#   test/comparison/  test_{scientist,coscientist}_params.{pdf,png,svg}
#
# Usage: source("NRC_heatmap_analysis.R") in RStudio, or
#        Rscript NRC_heatmap_analysis.R   (from the R/ directory)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ape)
  library(pheatmap)
  library(svglite)
})

################################################################################
# CONFIGURATION
################################################################################

ALL_PARAMS <- c("ipTM_LCB", "Sum_CONTACTS", "BSA_INT_PROTO", "SD_THETA_ROT",
                "S_PROTO", "D_APEX", "THETA_APEX_mean", "L_APEX_mean",
                "H_ABS_mean", "MU_H_mean", "D_MHD_P_mean")

SCIENTIST_ONLY      <- c("ipTM_LCB", "Sum_CONTACTS", "S_PROTO", "L_APEX_mean")
COSCIENTIST_ONLY    <- c("D_MHD_P_mean", "THETA_APEX_mean", "H_ABS_mean", "MU_H_mean")
OVERLAPPING         <- c("BSA_INT_PROTO", "SD_THETA_ROT", "D_APEX")

# Y-axis ordering: Scientist -> Both -> Co-Scientist
PARAM_ORDER <- c(SCIENTIST_ONLY, OVERLAPPING, COSCIENTIST_ONLY)
PARAM_SOURCE <- setNames(
  c(rep("Scientist",    length(SCIENTIST_ONLY)),
    rep("Both",         length(OVERLAPPING)),
    rep("Co-Scientist", length(COSCIENTIST_ONLY))),
  PARAM_ORDER
)
SOURCE_ORDER <- c("Scientist", "Both", "Co-Scientist")
SOURCE_COLOURS <- c(
  "Scientist"    = "#264653",
  "Both"         = "#2a9d8f",
  "Co-Scientist" = "#e76f51"
)

# For comparison folder: include overlapping in BOTH
SCIENTIST_PLUS_OVERLAP   <- c(SCIENTIST_ONLY, OVERLAPPING)
COSCIENTIST_PLUS_OVERLAP <- c(COSCIENTIST_ONLY, OVERLAPPING)

# Metric categories (for annotation_row colour bar)
metric_categories <- tibble(
  metric = ALL_PARAMS,
  category = case_when(
    metric %in% c("ipTM_LCB", "Sum_CONTACTS", "BSA_INT_PROTO") ~ "Confidence & Interface",
    metric %in% c("SD_THETA_ROT", "S_PROTO", "D_APEX")         ~ "Ring Symmetry",
    metric %in% c("THETA_APEX_mean", "L_APEX_mean",
                   "H_ABS_mean", "MU_H_mean")                   ~ "CC-Domain Helix",
    metric == "D_MHD_P_mean"                                     ~ "Motif Distances",
    TRUE ~ "Other"
  )
)

category_colours <- c(
  "Confidence & Interface" = "#003049",
  "Ring Symmetry"          = "#d52a28",
  "CC-Domain Helix"        = "#f78e1e",
  "Motif Distances"        = "#ffd641"
)

HEATMAP_COLORS <- c("#1C39BB", "white", "#ED2939")
REMOVE_FROM_TEST <- c("nrc0", "slnrc0_sa")

################################################################################
# HELPER FUNCTIONS
################################################################################

read_clade_trees <- function(tree_dir) {
  tree_files <- list.files(tree_dir, pattern = "\\.tree$", full.names = TRUE)
  clade_map <- map_dfr(tree_files, function(f) {
    clade_name <- tools::file_path_sans_ext(basename(f))
    tree <- read.tree(f)
    tips <- tree$tip.label
    tibble(ID_original = tips, clade = clade_name)
  })
  return(clade_map)
}

# Format raw values per metric (matching old script formatting)
fmt_raw <- function(raw_mat) {
  fmt <- matrix(sprintf("%.1f", raw_mat), nrow = nrow(raw_mat), dimnames = dimnames(raw_mat))
  if ("H_ABS_mean" %in% rownames(fmt))
    fmt["H_ABS_mean", ] <- sprintf("%.2f", raw_mat["H_ABS_mean", ])
  if ("BSA_INT_PROTO" %in% rownames(fmt))
    fmt["BSA_INT_PROTO", ] <- sprintf("%.0f", raw_mat["BSA_INT_PROTO", ])
  if ("S_PROTO" %in% rownames(fmt))
    fmt["S_PROTO", ] <- sprintf("%.2f", raw_mat["S_PROTO", ])
  if ("MU_H_mean" %in% rownames(fmt))
    fmt["MU_H_mean", ] <- sprintf("%.2f", raw_mat["MU_H_mean", ])
  if ("ipTM_LCB" %in% rownames(fmt))
    fmt["ipTM_LCB", ] <- sprintf("%.2f", raw_mat["ipTM_LCB", ])
  return(fmt)
}

# Save a pheatmap to PDF, PNG, and SVG
save_pheatmap <- function(p, filepath_base, w = 15, h = 15) {
  ggsave(paste0(filepath_base, ".pdf"), p, width = w, height = h)
  ggsave(paste0(filepath_base, ".png"), p, width = w, height = h, dpi = 300)
  tryCatch(
    ggsave(paste0(filepath_base, ".svg"), p, width = w, height = h, device = svglite),
    error = function(e) cat("SVG save skipped:", e$message, "\n")
  )
}

# Order a (metrics x clades) matrix by PARAM_ORDER and return the order +
# gap positions used for pheatmap's gaps_row (between Scientist/Both/Co-Scientist).
order_metrics_by_source <- function(mat_t) {
  metric_order_use <- intersect(PARAM_ORDER, rownames(mat_t))
  mat_t <- mat_t[metric_order_use, , drop = FALSE]

  source_vec <- PARAM_SOURCE[metric_order_use]
  # gaps go after positions where source changes
  gaps <- which(diff(as.integer(factor(source_vec, levels = SOURCE_ORDER))) != 0)

  list(
    mat = mat_t,
    order = metric_order_use,
    gaps = gaps,
    source = source_vec
  )
}

# Build the row annotation data frame (Source + structural Category) for pheatmap
build_row_annotation <- function(metric_order_use) {
  data.frame(
    Source = factor(PARAM_SOURCE[metric_order_use], levels = SOURCE_ORDER),
    Category = metric_categories$category[match(metric_order_use, metric_categories$metric)],
    row.names = metric_order_use,
    check.names = FALSE
  )
}

################################################################################
# MAIN ANALYSIS FUNCTION
################################################################################

analyze_dataset <- function(data_path, tree_dir, output_base,
                            dataset_name, remove_samples = NULL,
                            individual_mode = FALSE) {

  cat("\n", strrep("=", 70), "\n")
  cat("Analyzing dataset:", toupper(dataset_name),
      if (individual_mode) "(individual mode)" else "(clade mode)", "\n")
  cat(strrep("=", 70), "\n\n")

  # ---- 1. Read data ----
  raw <- read_excel(data_path)
  cat("Loaded", nrow(raw), "samples x", ncol(raw), "columns\n")

  # Remove specified samples
  if (!is.null(remove_samples)) {
    before <- nrow(raw)
    raw <- raw %>% filter(!tolower(ID) %in% tolower(remove_samples))
    cat("Removed", before - nrow(raw), "samples:", paste(remove_samples, collapse = ", "), "\n")
    cat("Remaining:", nrow(raw), "samples\n")
  }

  # ---- 2. Assign clade membership ----
  if (individual_mode) {
    # Each row is its own "clade" (used for the test/reference set)
    raw <- raw %>% mutate(clade = ID)
    cat("Individual mode: using each entry as its own group (", nrow(raw), "entries)\n")
  } else {
    # Read clade trees and join
    clade_map <- read_clade_trees(tree_dir)

    # Match by normalised ID (lowercase, dots -> underscores)
    raw <- raw %>%
      mutate(ID_norm = tolower(gsub("\\.", "_", ID)))
    clade_map <- clade_map %>%
      mutate(ID_norm = tolower(gsub("\\.", "_", ID_original)))

    raw <- raw %>%
      left_join(clade_map %>% select(ID_norm, clade), by = "ID_norm") %>%
      filter(!is.na(clade))

    cat("Matched to clades:", nrow(raw), "samples\n")
    cat("Clade distribution:\n")
    print(table(raw$clade))
  }

  # ---- 3. Available parameters ----
  score_cols <- intersect(ALL_PARAMS, colnames(raw))
  # Drop params that are all NA
  score_cols <- score_cols[sapply(score_cols, function(p) sum(!is.na(raw[[p]])) > 0)]
  cat("Available parameters:", paste(score_cols, collapse = ", "), "\n")

  # ---- 4. Compute clade means ----
  # Clamp negatives to 0
  raw_clamped <- raw %>%
    mutate(across(all_of(score_cols), ~ pmax(.x, 0)))

  clade_means <- raw_clamped %>%
    group_by(clade) %>%
    summarise(
      n_members = n(),
      across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  # Prepare matrix
  clade_means_mat <- clade_means %>%
    select(-n_members) %>%
    column_to_rownames("clade") %>%
    as.matrix()

  # Drop zero-variance columns
  clade_means_mat <- clade_means_mat[, apply(clade_means_mat, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  clade_means_scaled <- scale(clade_means_mat)  # Z-score columns

  # Annotation: n_members per clade
  annotation_row_clades <- clade_means %>%
    select(clade, n_members) %>%
    column_to_rownames("clade")

  # Annotation: metric categories
  annotation_col_metrics <- metric_categories %>%
    filter(metric %in% colnames(clade_means_scaled)) %>%
    column_to_rownames("metric")

  ann_colors <- list(
    Category = category_colours,
    category = category_colours,   # kept for back-compat with old annotation var
    Source   = SOURCE_COLOURS,
    n_members = c("white", "#313131")
  )

  # Create output directories
  heatmap_dir    <- file.path(output_base, "heatmaps")
  comparison_dir <- file.path(output_base, "comparison")
  for (d in c(heatmap_dir, comparison_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }

  # ==========================================================================
  # A. HEATMAPS - pheatmap, metrics manually ordered Scientist->Both->Co-Scientist
  # ==========================================================================
  cat("\n[A] Generating heatmaps...\n")

  # Transposed matrix: metrics x clades
  mat_t_all_full  <- t(clade_means_scaled)
  raw_all_t_full  <- t(clade_means_mat[, colnames(clade_means_scaled), drop = FALSE])

  # Manual row order by source (Scientist -> Both -> Co-Scientist)
  ord_all     <- order_metrics_by_source(mat_t_all_full)
  mat_t_all   <- ord_all$mat
  raw_all_t   <- raw_all_t_full[ord_all$order, , drop = FALSE]
  raw_all_fmt <- fmt_raw(raw_all_t)
  gaps_all    <- ord_all$gaps
  ann_row_all <- build_row_annotation(ord_all$order)

  # A1. All parameters - Z-scored with raw annotations
  p_heatmap_zraw <- pheatmap(
    mat_t_all,
    main = paste0(toupper(dataset_name), " Set \u2014 Clade Mean Scores (Z-scored, raw annotations)"),
    cluster_rows = FALSE,                         # manual order
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    gaps_row = gaps_all,                          # gaps between Scientist/Both/Co-Sci
    annotation_row = ann_row_all,
    annotation_col = annotation_row_clades,
    annotation_colors = ann_colors,
    color = colorRampPalette(HEATMAP_COLORS)(100),
    breaks = seq(-3, 3, length.out = 101),
    fontsize = 11,
    fontsize_row = 11,
    fontsize_col = 11,
    angle_col = 45,
    border_color = "white",
    display_numbers = raw_all_fmt,
    number_color = "black",
    cellwidth = 30, cellheight = 30
  )
  save_pheatmap(p_heatmap_zraw, file.path(heatmap_dir, paste0(dataset_name, "_clade_heatmap_zraw")))
  cat("  Saved all-parameters heatmap\n")

  # A2. Z-scored only (no numbers)
  p_heatmap_z <- pheatmap(
    mat_t_all,
    main = paste0(toupper(dataset_name), " Set \u2014 Clade Mean Scores (Z-scored)"),
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    clustering_method = "ward.D2",
    gaps_row = gaps_all,
    annotation_row = ann_row_all,
    annotation_col = annotation_row_clades,
    annotation_colors = ann_colors,
    color = colorRampPalette(HEATMAP_COLORS)(100),
    breaks = seq(-3, 3, length.out = 101),
    fontsize = 10,
    fontsize_row = 11,
    fontsize_col = 11,
    angle_col = 45,
    border_color = "white",
    display_numbers = FALSE,
    cellwidth = 25, cellheight = 25
  )
  save_pheatmap(p_heatmap_z, file.path(heatmap_dir, paste0(dataset_name, "_clade_heatmap")))

  # ==========================================================================
  # B. COMPARISON HEATMAPS (overlapping params in BOTH Scientist & Co-Scientist)
  # ==========================================================================
  cat("\n[B] Generating comparison heatmaps...\n")

  comparison_sets <- list(
    "Scientist + Overlapping"    = intersect(SCIENTIST_PLUS_OVERLAP, colnames(clade_means_scaled)),
    "Co-Scientist + Overlapping" = intersect(COSCIENTIST_PLUS_OVERLAP, colnames(clade_means_scaled))
  )

  for (set_name in names(comparison_sets)) {
    cols <- comparison_sets[[set_name]]
    if (length(cols) == 0) next

    sub_scaled <- clade_means_scaled[, cols, drop = FALSE]
    sub_raw <- clade_means_mat[, cols, drop = FALSE]

    # Manual row order by source
    ord_sub     <- order_metrics_by_source(t(sub_scaled))
    sub_mat_t   <- ord_sub$mat
    sub_raw_t   <- t(sub_raw)[ord_sub$order, , drop = FALSE]
    sub_raw_fmt <- fmt_raw(sub_raw_t)
    sub_ann_row <- build_row_annotation(ord_sub$order)

    tag <- if (grepl("Co-Sci", set_name)) "coscientist" else "scientist"

    p_comp <- pheatmap(
      sub_mat_t,
      main = paste0(toupper(dataset_name), " Set \u2014 ", set_name, " Metrics (Z-scored)"),
      cluster_rows = FALSE,
      cluster_cols = TRUE,
      clustering_method = "ward.D2",
      gaps_row = ord_sub$gaps,
      annotation_row = sub_ann_row,
      annotation_col = annotation_row_clades,
      annotation_colors = ann_colors,
      color = colorRampPalette(HEATMAP_COLORS)(100),
      breaks = seq(-3, 3, length.out = 101),
      fontsize = 11,
      fontsize_row = 11,
      fontsize_col = 11,
      angle_col = 45,
      border_color = "white",
      display_numbers = sub_raw_fmt,
      number_color = "black",
      cellwidth = 30, cellheight = 30
    )
    save_pheatmap(p_comp, file.path(comparison_dir, paste0(dataset_name, "_", tag, "_params")))
    cat("  Saved", set_name, "heatmap\n")
  }

  cat("\n", strrep("=", 70), "\n")
  cat("Analysis complete for:", toupper(dataset_name), "\n")
  cat(strrep("=", 70), "\n\n")

  return(list(
    raw = raw,
    clade_means = clade_means,
    clade_means_scaled = clade_means_scaled
  ))
}

################################################################################
# EXECUTE ANALYSIS
################################################################################

# Try to set working directory to script location
tryCatch({
  script_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  setwd(script_dir)
}, error = function(e) {
  cat("Not in RStudio; using current working directory.\n")
  cat("Make sure you run this from the R/ directory.\n")
})

cat("Working directory:", getwd(), "\n\n")

# Data paths (relative to R/)
main_path <- "main/resistosome_analysis_summary.xlsx"
test_path <- "test/resistosome_analysis_summary.xlsx"
tree_dir  <- "../phylo/clades/"

# ---- Main set ----
if (file.exists(main_path)) {
  main_results <- analyze_dataset(
    data_path    = main_path,
    tree_dir     = tree_dir,
    output_base  = "main/",
    dataset_name = "main"
  )
} else {
  cat("Main data not found at:", normalizePath(main_path, mustWork = FALSE), "\n")
}

# ---- Test set (remove NRC0 and SlNRC0-Sa) ----
if (file.exists(test_path)) {
  tryCatch({
    test_results <- analyze_dataset(
      data_path       = test_path,
      tree_dir        = tree_dir,
      output_base     = "test/",
      dataset_name    = "test",
      remove_samples  = REMOVE_FROM_TEST,
      individual_mode = TRUE
    )
  }, error = function(e) {
    cat("Test set analysis error:", e$message, "\n")
  })
} else {
  cat("Test data not found at:", normalizePath(test_path, mustWork = FALSE), "\n")
}

cat("\n\nAll analyses complete!\n")
