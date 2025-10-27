#!/usr/bin/env Rscript

###############################################################################
# H2O AutoML training pipeline driven by config.yml.
#
# Usage:
#   Rscript scripts/03_h2o_train.R [dataset]
#
# The script:
#   * Loads environmental rasters for the dataset.
#   * Applies deterministic Kendall τ pruning (cached in plans/).
#   * Runs H2O AutoML on each occurrence file (presence/background balance).
#   * Performs spatial block cross-validation and exports metrics.
#   * Generates raster predictions and auxiliary artifacts.
###############################################################################

suppressPackageStartupMessages({
  library(terra)
  library(h2o)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(sf)
  library(blockCV)
  library(tools)
  library(fs)
})

source("R/utils_io.R")
source("R/utils_repro.R")
source("R/utils_kendall.R")
source("R/utils_h2o.R")
source("R/utils_plot.R")

cfg <- get_config()
args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else cfg$runtime$dataset
stopifnot(dataset %in% names(cfg$datasets))

dataset_cfg <- cfg$datasets[[dataset]]
occ_files <- vapply(dataset_cfg$occ_files, function(p) resolve_path(cfg, p), character(1))
env_dir   <- resolve_path(cfg, dataset_cfg$env_dir)
plan_dir  <- resolve_path(cfg, cfg$paths$plans, dataset)
results_root <- resolve_path(cfg, dataset_cfg$results$h2o)
ensure_dir(plan_dir)
ensure_dir(results_root)

with_repro_context({
  env_files <- fs::dir_ls(env_dir, glob = "*.tif")
  if (length(env_files) == 0) {
    stop(sprintf("No raster files (*.tif) found in %s", env_dir))
  }
  env_names <- make.names(fs::path_ext_remove(fs::path_file(env_files)), unique = TRUE)
  env_stack <- terra::rast(env_files)
  names(env_stack) <- env_names

  kendall_keep_path <- fs::path(plan_dir, "keepvars.csv")
  kendall_drop_path <- fs::path(plan_dir, "kendall_dropped_vars.csv")
  kendall_heatmap_path <- fs::path(plan_dir, "kendall_heatmap.png")
  kendall_matrix_path <- fs::path(plan_dir, "kendall_matrix.rds")

  if (fs::file_exists(kendall_keep_path)) {
    keep_vars <- readr::read_csv(kendall_keep_path, show_col_types = FALSE)$variable
    drop_vars <- if (fs::file_exists(kendall_drop_path)) {
      readr::read_csv(kendall_drop_path, show_col_types = FALSE)$variable
    } else character(0)
    kt <- if (fs::file_exists(kendall_matrix_path)) readRDS(kendall_matrix_path) else NULL
  } else {
    message("Computing Kendall correlation matrix ...")
    kt <- compute_kendall_matrix(env_stack, cfg$modeling$kendall_sample)
    drop_vars <- prune_by_kendall(kt, cutoff = cfg$modeling$kendall_cutoff)
    keep_vars <- setdiff(colnames(kt), drop_vars)
    readr::write_csv(tibble::tibble(variable = keep_vars), kendall_keep_path)
    readr::write_csv(tibble::tibble(variable = drop_vars), kendall_drop_path)
    saveRDS(kt, kendall_matrix_path)
    plot_kendall_heatmap(kt[keep_vars, keep_vars, drop = FALSE],
                         cfg$modeling$kendall_cutoff,
                         kendall_heatmap_path)
  }

  if (length(keep_vars) == 0) {
    stop("Kendall screening removed all variables; check cutoff in config.yml")
  }

  env_stack <- env_stack[[keep_vars]]
  message(sprintf("Using %d environmental variables after Kendall pruning.", length(keep_vars)))

  h2o::h2o.init()

  for (occ in occ_files) {
    if (!fs::file_exists(occ)) {
      warning(sprintf("Occurrence file missing: %s", occ))
      next
    }
    sp <- read_csv_safely(occ)
    stopifnot(all(c("lon", "lat") %in% colnames(sp)))

    dset <- fs::path_ext_remove(fs::path_file(occ))
    species_dir <- fs::path(results_root, dset)
    ensure_dir(species_dir)

    message("====================================================")
    message(sprintf("Dataset: %s", dset))

    pres <- sp %>% dplyr::select(lon, lat) %>% dplyr::mutate(pa = 1)
    bg_pts <- terra::spatSample(env_stack[[1]], size = nrow(pres),
                                method = "random", as.points = TRUE)
    bg_coords <- terra::crds(bg_pts)
    bg_df <- tibble::tibble(lon = bg_coords[, 1], lat = bg_coords[, 2], pa = 0)

    df_sp <- dplyr::bind_rows(pres, bg_df)
    pts <- terra::vect(df_sp, geom = c("lon", "lat"), crs = terra::crs(env_stack))
    vals <- terra::extract(env_stack, pts)[, -1, drop = FALSE]
    df <- dplyr::bind_cols(df_sp, as.data.frame(vals)) %>% tidyr::drop_na()

    sf_pts <- sf::st_as_sf(df, coords = c("lon", "lat"), crs = terra::crs(env_stack))

    folds_obj <- NULL
    if ("cv_spatial" %in% getNamespaceExports("blockCV")) {
      folds_obj <- blockCV::cv_spatial(x = sf_pts, column = "pa",
                                       size = cfg$modeling$block_km * 1000,
                                       k = cfg$modeling$n_folds)
      if (!is.null(folds_obj$folds_ids)) {
        df$fold_id <- as.integer(folds_obj$folds_ids)
      } else if (!is.null(folds_obj$folds_list)) {
        fold_id <- rep(NA_integer_, nrow(df))
        for (k in seq_along(folds_obj$folds_list)) {
          fold_id[folds_obj$folds_list[[k]]] <- k
        }
        df$fold_id <- fold_id
      } else {
        stop("cv_spatial(): cannot find fold ids; check blockCV version.")
      }
    } else {
      sb <- blockCV::spatialBlock(
        speciesData = sf_pts, species = "pa",
        theRange = cfg$modeling$block_km * 1000,
        k = cfg$modeling$n_folds,
        selection = "systematic",
        biomod2Format = FALSE,
        showBlocks = FALSE
      )
      fold_id <- rep(NA_integer_, nrow(df))
      for (k in seq_along(sb$folds)) {
        fold_id[sb$folds[[k]]$test] <- k
      }
      df$fold_id <- fold_id
    }
    stopifnot(!any(is.na(df$fold_id)))

    hf <- h2o::as.h2o(df)
    hf[["pa"]] <- h2o::asfactor(hf[["pa"]])
    aml <- h2o::h2o.automl(
      x = keep_vars,
      y = "pa",
      training_frame = hf,
      max_models = cfg$modeling$automl_models,
      seed = cfg$project$seed
    )
    leader <- aml@leader
    model_path <- h2o::h2o.saveModel(leader, path = species_dir, force = TRUE)
    message(sprintf("Saved leader model: %s", model_path))

    perf <- h2o::h2o.performance(leader, newdata = hf)
    metrics <- tibble::tibble(
      dataset = dset,
      rmse = h2o::h2o.rmse(perf),
      auc = h2o::h2o.auc(perf),
      logloss = h2o::h2o.logloss(perf)
    )
    write_csv_safely(metrics, fs::path(species_dir, sprintf("metrics_in_sample_%s.csv", dset)))

    varimp <- tryCatch(as.data.frame(h2o::h2o.varimp(leader)), error = function(e) NULL)
    if (!is.null(varimp) && nrow(varimp) > 0) {
      varimp$dataset <- dset
      write_csv_safely(varimp, fs::path(species_dir, sprintf("varimp_%s.csv", dset)))
      pdf(fs::path(species_dir, sprintf("varimp_%s.pdf", dset)))
      print(h2o::h2o.varimp_plot(leader, num_of_features = min(10, nrow(varimp))))
      dev.off()
      pp_vars <- head(keep_vars, cfg$modeling$partial_vars)
      suppressWarnings({
        pdf(fs::path(species_dir, sprintf("partial_%s.pdf", dset)))
        print(h2o::h2o.partialPlot(object = leader, data = hf, cols = pp_vars))
        dev.off()
      })
    }

    cv_rows <- list()
    for (k in sort(unique(df$fold_id))) {
      trn <- df[df$fold_id != k, , drop = FALSE]
      tst <- df[df$fold_id == k, , drop = FALSE]
      if (nrow(trn) < 50 || nrow(tst) < 20) {
        next
      }
      hf_trn <- h2o::as.h2o(trn)
      hf_trn[["pa"]] <- h2o::asfactor(hf_trn[["pa"]])
      hf_tst <- h2o::as.h2o(tst)
      hf_tst[["pa"]] <- h2o::asfactor(hf_tst[["pa"]])

      aml_cv <- h2o::h2o.automl(
        x = keep_vars,
        y = "pa",
        training_frame = hf_trn,
        max_models = cfg$modeling$cv_automl_models,
        seed = cfg$project$seed + k
      )
      leader_cv <- aml_cv@leader
      perf_cv <- h2o::h2o.performance(leader_cv, newdata = hf_tst)

      cv_rows[[length(cv_rows) + 1]] <- tibble::tibble(
        dataset = dset,
        fold = k,
        auc = h2o::h2o.auc(perf_cv),
        rmse = h2o::h2o.rmse(perf_cv),
        logloss = h2o::h2o.logloss(perf_cv)
      )
    }

    if (length(cv_rows) > 0) {
      cv_df <- dplyr::bind_rows(cv_rows)
      cv_sum <- cv_df %>%
        summarise(
          dataset = dset,
          folds = dplyr::n(),
          auc_mean = mean(auc, na.rm = TRUE),
          auc_sd = sd(auc, na.rm = TRUE),
          rmse_mean = mean(rmse, na.rm = TRUE),
          rmse_sd = sd(rmse, na.rm = TRUE),
          logloss_mean = mean(logloss, na.rm = TRUE),
          logloss_sd = sd(logloss, na.rm = TRUE)
        )
      write_csv_safely(cv_df, fs::path(species_dir, sprintf("spatialCV_folds_%s.csv", dset)))
      write_csv_safely(cv_sum, fs::path(species_dir, sprintf("spatialCV_summary_%s.csv", dset)))

      g1 <- ggplot(cv_df, aes(x = factor(fold), y = auc)) +
        geom_col() +
        labs(x = "Fold", y = "AUC (spatial CV)", title = sprintf("Spatial CV — %s", dset)) +
        theme_minimal(base_size = 11)
      ggplot2::ggsave(fs::path(species_dir, sprintf("spatialCV_auc_%s.png", dset)),
                      g1, width = 5.5, height = 3.6, dpi = 300)
    }

    pred_tif <- fs::path(species_dir, sprintf("prediction_%s.tif", dset))
    predict_raster_h2o(env_stack, leader, pred_tif, batch = cfg$modeling$prediction_batch)
    message(sprintf("Saved prediction raster: %s", pred_tif))
  }

  try(h2o::h2o.shutdown(prompt = FALSE), silent = TRUE)
}, seed = cfg$project$seed, single_core = cfg$project$single_core)

message("Done.")

