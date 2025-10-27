#!/usr/bin/env Rscript

###############################################################################
# Compare H2O AutoML predictions against SSDM outputs using config.yml paths.
# Generates per-dataset raster comparisons, metrics, and figures.
###############################################################################

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(fs)
})

source("R/utils_io.R")
source("R/utils_repro.R")

cfg <- get_config()

ssdm_dir <- resolve_path(cfg, cfg$paths$results$ssdm)
h2o_dir  <- resolve_path(cfg, cfg$paths$results$h2o)
out_dir  <- resolve_path(cfg, cfg$paths$results$compare)
ensure_dir(out_dir)
ensure_dir(fs::path(out_dir, "01_between_methods", "rasters"))
ensure_dir(fs::path(out_dir, "01_between_methods", "plots"))
ensure_dir(fs::path(out_dir, "02_temporal", "rasters"))
ensure_dir(fs::path(out_dir, "02_temporal", "plots"))
ensure_dir(fs::path(out_dir, "03_maps"))
ensure_dir(fs::path(out_dir, "04_panels"))

with_repro_context({
  cat(
    "Folders:\n",
    " 01_between_methods: pixel-wise SSDM vs H2O (diff rasters + metrics plots)\n",
    " 02_temporal: Δ (A−B) rasters, gain/stable/loss rasters\n",
    " 03_maps: base suitability maps per dataset & method\n",
    " 04_panels: combined panels per dataset\n",
    file = fs::path(out_dir, "README.txt")
  )

  parse_meta <- function(fp, method_hint) {
    nm <- fs::path_file(fp)
    dataset <- stringr::str_match(nm, "(E[0-9]+[AB])")[, 2]
    tibble::tibble(file = fp, dataset = dataset, method = method_hint) %>%
      tidyr::drop_na(dataset)
  }

  ssdm_files <- fs::dir_ls(ssdm_dir, glob = "**/*.tif")
  h2o_files  <- fs::dir_ls(h2o_dir, glob = "**/prediction_*.tif")

  meta_ssdm <- dplyr::bind_rows(lapply(ssdm_files, parse_meta, method_hint = "SSDM"))
  meta_h2o  <- dplyr::bind_rows(lapply(h2o_files,  parse_meta, method_hint = "H2O"))

  meta_all <- dplyr::full_join(meta_ssdm, meta_h2o, by = "dataset",
                               suffix = c("_ssdm", "_h2o"))
  datasets_all <- stats::na.omit(meta_all$dataset)
  if (length(datasets_all) == 0) {
    stop("No paired datasets found between SSDM and H2O rasters.")
  }

  align_to <- function(r1, r2, categorical = FALSE) {
    if (!terra::compareGeom(r1, r2, stopOnError = FALSE)) {
      r2 <- terra::project(r2, terra::crs(r1),
                           method = if (categorical) "near" else "bilinear")
      r2 <- terra::resample(r2, r1, method = if (categorical) "near" else "bilinear")
      r2 <- terra::crop(r2, r1)
    }
    m <- !is.na(r1) & !is.na(r2)
    r1 <- terra::mask(r1, m, maskvalues = 0)
    r2 <- terra::mask(r2, m, maskvalues = 0)
    list(r1 = r1, r2 = r2)
  }

  r_to_df_full <- function(r) {
    as.data.frame(r, xy = TRUE, na.rm = TRUE) |>
      stats::setNames(c("x", "y", "val"))
  }

  plot_raster_continuous <- function(r, title, out_png, center0 = FALSE) {
    df <- r_to_df_full(r)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = val)) +
      ggplot2::geom_raster() +
      ggplot2::coord_equal() +
      ggplot2::labs(title = title, x = NULL, y = NULL, fill = NULL) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0))
    if (center0) {
      lim <- max(abs(range(df$val, na.rm = TRUE)))
      p <- p + ggplot2::scale_fill_gradientn(
        colors = c("#2c7fb8", "#ffffbf", "#d7191c"),
        limits = c(-lim, lim), na.value = NA
      )
    } else {
      p <- p + ggplot2::scale_fill_viridis_c(na.value = NA)
    }
    ggplot2::ggsave(out_png, p, width = 8, height = 6, dpi = 300)
  }

  plot_raster_discrete <- function(r, title, out_png) {
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(df) <- c("x", "y", "val")
    df$val <- factor(df$val, levels = c(-1, 0, 1),
                     labels = c("Loss", "Stable", "Gain"))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = val)) +
      ggplot2::geom_raster() +
      ggplot2::coord_equal() +
      ggplot2::scale_fill_manual(
        values = c("Loss" = "#b2182b", "Stable" = "#f7f7f7", "Gain" = "#2166ac"),
        drop = FALSE
      ) +
      ggplot2::labs(title = title, x = NULL, y = NULL, fill = NULL) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(hjust = 0))
    ggplot2::ggsave(out_png, p, width = 8, height = 6, dpi = 300)
  }

  jaccard_binary <- function(b1, b2) {
    inter <- terra::global(b1 & b2, "sum", na.rm = TRUE)[[1]]
    union <- terra::global(b1 | b2, "sum", na.rm = TRUE)[[1]]
    ifelse(union == 0, NA_real_, inter / union)
  }

  quick_bar <- function(df, metric, ylab, outpng) {
    df$ord <- factor(df$dataset, levels = sort(unique(df$dataset), decreasing = TRUE))
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = "ord", y = metric)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = ylab,
                    title = paste("Between-method:", ylab)) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0))
    ggplot2::ggsave(outpng, p, width = 7.5, height = 5, dpi = 300)
  }

  bm_metrics <- list()

  for (ds in datasets_all) {
    f_ssdm <- meta_all %>% dplyr::filter(dataset == ds) %>% dplyr::pull(file_ssdm)
    f_h2o  <- meta_all %>% dplyr::filter(dataset == ds) %>% dplyr::pull(file_h2o)
    if (length(f_ssdm) * length(f_h2o) == 0 || anyNA(c(f_ssdm, f_h2o))) {
      message(sprintf("Skipping %s (missing pair)", ds))
      next
    }
    message(sprintf("Between-method comparison: %s", ds))

    r_ssdm <- terra::rast(f_ssdm)
    r_h2o  <- terra::rast(f_h2o)
    al <- align_to(r_ssdm, r_h2o, categorical = FALSE)
    r_ssdm <- al$r1
    r_h2o  <- al$r2

    plot_raster_continuous(r_ssdm, sprintf("%s — SSDM suitability", ds),
                           fs::path(out_dir, "03_maps", sprintf("%s_SSDM.png", ds)))
    plot_raster_continuous(r_h2o, sprintf("%s — H2O suitability", ds),
                           fs::path(out_dir, "03_maps", sprintf("%s_H2O.png", ds)))

    v1 <- terra::values(r_ssdm, mat = FALSE)
    v2 <- terra::values(r_h2o, mat = FALSE)
    diff_r <- r_h2o - r_ssdm
    plot_raster_continuous(diff_r, sprintf("%s — H2O minus SSDM", ds),
                           fs::path(out_dir, "01_between_methods", "rasters",
                                    sprintf("%s_diff.png", ds)), center0 = TRUE)

    df_metrics <- tibble::tibble(
      dataset = ds,
      cor = stats::cor(v1, v2, use = "complete.obs"),
      mae = mean(abs(v1 - v2), na.rm = TRUE),
      rmse = sqrt(mean((v1 - v2)^2, na.rm = TRUE))
    )
    bm_metrics[[length(bm_metrics) + 1]] <- df_metrics

  }

  if (length(bm_metrics) > 0) {
    bm_df <- dplyr::bind_rows(bm_metrics)
    write_csv_safely(bm_df, fs::path(out_dir, "01_between_methods", "metrics.csv"))
    quick_bar(bm_df, "cor", "Pearson correlation", fs::path(out_dir, "01_between_methods", "plots", "correlation.png"))
    quick_bar(bm_df, "mae", "MAE (suitability)", fs::path(out_dir, "01_between_methods", "plots", "mae.png"))
    quick_bar(bm_df, "rmse", "RMSE (suitability)", fs::path(out_dir, "01_between_methods", "plots", "rmse.png"))
  }
}, seed = cfg$project$seed, single_core = cfg$project$single_core)

