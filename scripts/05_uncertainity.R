#!/usr/bin/env Rscript

###############################################################################
# Variable importance comparison (before vs after) using H2O AutoML outputs.
# Reads species definitions from config.yml under `uncertainty$species`.
###############################################################################

suppressPackageStartupMessages({
  library(h2o)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(patchwork)
  library(readr)
  library(fs)
})

source("R/utils_io.R")
source("R/utils_repro.R")

cfg <- get_config()

%||% <- function(x, y) if (!is.null(x)) x else y

if (is.null(cfg$uncertainty$species) || length(cfg$uncertainty$species) == 0) {
  stop("No species configured under `uncertainty$species` in config.yml")
}

species_dirs <- lapply(cfg$uncertainty$species, function(paths) {
  list(
    before = resolve_path(cfg, paths$before),
    after  = resolve_path(cfg, paths$after)
  )
})
names(species_dirs) <- names(cfg$uncertainty$species)

out_dir <- resolve_path(cfg, cfg$uncertainty$output_dir)
ensure_dir(out_dir)

with_repro_context({
  h2o::h2o.init()

  BIO_LUT <- c(
    bio1  = "Annual Mean Temperature",
    bio2  = "Mean Diurnal Range",
    bio3  = "Isothermality",
    bio4  = "Temperature Seasonality",
    bio5  = "Max Temp of Warmest Month",
    bio6  = "Min Temp of Coldest Month",
    bio7  = "Temperature Annual Range",
    bio8  = "Mean Temp of Wettest Quarter",
    bio9  = "Mean Temp of Driest Quarter",
    bio10 = "Mean Temp of Warmest Quarter",
    bio11 = "Mean Temp of Coldest Quarter",
    bio12 = "Annual Precipitation",
    bio13 = "Precipitation of Wettest Month",
    bio14 = "Precipitation of Driest Month",
    bio15 = "Precipitation Seasonality",
    bio16 = "Precipitation of Wettest Quarter",
    bio17 = "Precipitation of Driest Quarter",
    bio18 = "Precipitation of Warmest Quarter",
    bio19 = "Precipitation of Coldest Quarter"
  )

  pretty_var <- function(x) {
    x_chr <- as.character(x)
    lower <- tolower(x_chr)
    out <- x_chr
    map_idx <- match(lower, names(BIO_LUT))
    mapped <- BIO_LUT[map_idx]
    repl <- !is.na(map_idx)
    out[repl] <- mapped[repl]
    idx_year <- grepl("^(ndvi|evi)_?\\d{4}$", lower, perl = TRUE)
    if (any(idx_year)) {
      yr <- stringr::str_extract(x_chr[idx_year], "\\\d{4}$")
      base <- toupper(sub("_?\\\d{4}$", "", x_chr[idx_year], perl = TRUE))
      out[idx_year] <- paste0(base, " (", yr, ")")
    }
    out[grepl("^ndvi$", lower)] <- "NDVI"
    out[grepl("^evi$", lower)] <- "EVI"
    out[grepl("^elev|^elevation$", lower)] <- "Elevation (m)"
    out[grepl("^slope$", lower)] <- "Slope (°)"
    out[grepl("^aspect$", lower)] <- "Cos_Aspect (°)"
    out[grepl("^land\\s*cover$", lower)] <- "LandCover"
    out[grepl("^forest\\s*cover$", lower)] <- "ForestCover"
    still_raw <- out == x_chr
    if (any(still_raw)) {
      tmp <- gsub("_", " ", out[still_raw])
      tmp <- stringr::str_trim(tmp)
      tmp <- stringr::str_squish(tmp)
      tmp <- stringr::str_to_title(tmp)
      out[still_raw] <- tmp
    }
    out
  }

  pick_model <- function(dir_path) {
    if (!fs::dir_exists(dir_path)) {
      stop(sprintf("Directory not found: %s", dir_path))
    }
    models <- fs::dir_ls(dir_path, type = "directory")
    stacked <- models[grepl("StackedEnsemble", fs::path_file(models))]
    if (length(stacked)) {
      return(stacked[[1]])
    }
    fallback <- models[grepl("XGBoost|GBM|DRF|DeepLearning|GLM", fs::path_file(models))]
    if (length(fallback)) {
      return(fallback[[1]])
    }
    if (length(models)) {
      return(models[[1]])
    }
    NA_character_
  }

  combined_varimp <- function(model) {
    if (is.na(model@model_id)) {
      return(tibble::tibble(variable = character(), rel_imp = numeric()))
    }
    is_ens <- grepl("StackedEnsemble", model@algorithm, ignore.case = TRUE)
    if (!is_ens) {
      vi <- try(h2o::h2o.varimp(model), silent = TRUE)
      if (inherits(vi, "try-error") || is.null(vi)) {
        return(tibble::tibble(variable = character(), rel_imp = numeric()))
      }
      dplyr::as_tibble(vi) %>%
        dplyr::select(variable, relative_importance) %>%
        dplyr::group_by(variable) %>%
        dplyr::summarise(rel_imp = sum(relative_importance), .groups = "drop") %>%
        dplyr::mutate(rel_imp = 100 * rel_imp / sum(rel_imp))
    } else {
      base_ids <- model@model$base_models
      base_ids <- base_ids[!grepl("Metalearner", base_ids, ignore.case = TRUE)]
      parts <- lapply(base_ids, function(bid) {
        bm <- h2o::h2o.getModel(bid)
        vi <- try(h2o::h2o.varimp(bm), silent = TRUE)
        if (inherits(vi, "try-error") || is.null(vi)) {
          return(NULL)
        }
        dplyr::as_tibble(vi) %>%
          dplyr::select(variable, relative_importance)
      })
      parts <- Filter(Negate(is.null), parts)
      if (!length(parts)) {
        return(tibble::tibble(variable = character(), rel_imp = numeric()))
      }
      dplyr::bind_rows(parts) %>%
        dplyr::group_by(variable) %>%
        dplyr::summarise(rel_imp = sum(relative_importance), .groups = "drop") %>%
        dplyr::mutate(rel_imp = 100 * rel_imp / sum(rel_imp))
    }
  }

  plot_dumbbell <- function(wide_df, sp, top_k = Inf) {
    dd <- wide_df %>%
      dplyr::mutate(maxval = pmax(Before, After, na.rm = TRUE)) %>%
      dplyr::arrange(dplyr::desc(maxval))
    if (is.finite(top_k)) {
      dd <- dd %>% dplyr::slice_head(n = min(top_k, nrow(dd)))
    }
    dd$label <- pretty_var(dd$variable)
    dd$label <- factor(dd$label, levels = rev(dd$label))

    ggplot2::ggplot(dd) +
      ggplot2::geom_segment(ggplot2::aes(x = label, xend = label, y = Before, yend = After),
                            linewidth = 0.9, color = "grey65") +
      ggplot2::geom_point(ggplot2::aes(label, Before), size = 3.2, shape = 21,
                          fill = "#7A7A7A", color = "black", stroke = 0.25) +
      ggplot2::geom_point(ggplot2::aes(label, After), size = 3.2, shape = 21,
                          fill = "#1F77B4", color = "black", stroke = 0.25) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
      ggplot2::labs(title = paste0("Variable Importance — ", sp, " (Before vs After)"),
                    x = NULL, y = "Relative importance (%)") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 10.5)
      )
  }

  figure_height <- function(n_vars, base = 2.0, per_var = 0.33,
                            min_h = 6.0, max_h = 20.0) {
    h <- base + per_var * n_vars
    max(min_h, min(h, max_h))
  }

  all_wide <- list()
  all_plots <- list()

  for (sp in names(species_dirs)) {
    message(sprintf("Processing species %s", sp))
    dir_before <- species_dirs[[sp]][["before"]]
    dir_after  <- species_dirs[[sp]][["after"]]

    model_before_path <- pick_model(dir_before)
    model_after_path  <- pick_model(dir_after)

    if (anyNA(c(model_before_path, model_after_path))) {
      warning(sprintf("Skipping %s (missing model files)", sp))
      next
    }

    model_before <- h2o::h2o.loadModel(model_before_path)
    model_after  <- h2o::h2o.loadModel(model_after_path)

    vi_before <- combined_varimp(model_before)
    vi_after  <- combined_varimp(model_after)

    wide <- dplyr::full_join(vi_before, vi_after, by = "variable",
                             suffix = c("_before", "_after")) %>%
      dplyr::rename(Before = rel_imp_before, After = rel_imp_after) %>%
      tidyr::replace_na(list(Before = 0, After = 0))

    wide_with_labels <- wide %>% dplyr::mutate(label = pretty_var(variable))
    write_csv_safely(wide_with_labels, fs::path(out_dir, sprintf("varimp_%s.csv", sp)))

    plot <- plot_dumbbell(wide, sp, top_k = cfg$uncertainty$top_k %||% Inf)
    ggplot2::ggsave(fs::path(out_dir, sprintf("varimp_%s.png", sp)), plot,
                    width = 10, height = figure_height(nrow(wide)))

    all_wide[[sp]] <- wide
    all_plots[[sp]] <- plot
  }

  if (length(all_plots) > 0) {
    combined <- patchwork::wrap_plots(all_plots, ncol = 1)
    total_height <- sum(vapply(all_wide, function(df) figure_height(nrow(df)), numeric(1)))
    ggplot2::ggsave(fs::path(out_dir, "varimp_all_species.png"), combined,
                    width = 11, height = total_height)
  }

  try(h2o::h2o.shutdown(prompt = FALSE), silent = TRUE)
}, seed = cfg$project$seed, single_core = cfg$project$single_core)

