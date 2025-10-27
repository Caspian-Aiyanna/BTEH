#!/usr/bin/env Rscript

###############################################################################
# DBSCAN-based thinning for occurrence data with replicates.
#
# Usage:
#   Rscript scripts/02_dbscan_thin_degrees.R [dataset]
#
# Reads configuration from config.yml and writes thinned replicates to the
# dataset-specific occurrence directory defined therein.
###############################################################################

suppressPackageStartupMessages({
  library(dbscan)
  library(dplyr)
  library(readr)
  library(fs)
})

source("R/utils_io.R")
source("R/utils_repro.R")

cfg <- get_config()
args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else cfg$runtime$dataset
stopifnot(dataset %in% names(cfg$datasets))

dataset_cfg <- cfg$datasets[[dataset]]
input_dir <- resolve_path(cfg, cfg$paths$data[[cfg$thinning$input_subdir]], dataset)
output_dir <- resolve_path(cfg, cfg$paths$data[[cfg$thinning$output_subdir]], dataset)
ensure_dir(output_dir)

with_repro_context({
  csv_files <- fs::dir_ls(input_dir, glob = "*.csv")
  if (length(csv_files) == 0) {
    stop(sprintf("No CSV files found in %s", input_dir))
  }

  sample_with_dbscan <- function(df,
                                 eps = cfg$thinning$eps,
                                 minPts = cfg$thinning$min_pts,
                                 fraction = cfg$thinning$fraction,
                                 min_samples = cfg$thinning$min_samples) {
    db <- dbscan::dbscan(as.matrix(df[, c("lon", "lat")]), eps = eps, minPts = minPts)
    df$cluster_id <- db$cluster
    df %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::group_modify(~ {
        if (.y$cluster_id == 0) {
          dplyr::slice_sample(.x, prop = fraction)
        } else {
          n_to_sample <- max(min_samples, round(nrow(.x) * fraction))
          dplyr::slice_sample(.x, n = min(nrow(.x), n_to_sample))
        }
      }) %>%
      dplyr::ungroup()
  }

  for (csv in csv_files) {
    df <- readr::read_csv(csv, show_col_types = FALSE)
    stopifnot(all(c("lon", "lat") %in% colnames(df)))
    base <- fs::path_ext_remove(fs::path_file(csv))
    for (replicate in seq_len(cfg$thinning$replicates)) {
      set.seed(replicate)
      thinned <- sample_with_dbscan(df)
      out_file <- fs::path(output_dir, sprintf("%s_rep%d.csv", base, replicate))
      readr::write_csv(thinned, out_file)
      message(sprintf("Saved %s (%d points)", out_file, nrow(thinned)))
    }
  }
}, seed = cfg$project$seed, single_core = cfg$project$single_core)

