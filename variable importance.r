###############################################################################
# H2O Variable Importance — Before/After for E3, E4, E5 (ALL vars + FULL NAMES)
# - Finds leader (prefers StackedEnsemble*, else GBM/XGBoost/DRF/DeepLearning/GLM)
# - Combined VI across base learners (metalearner excluded), normalized to 100%
# - Plots ALL variables with FULL names; dynamic height; writes CSVs with labels
#
# Requires: h2o, dplyr, ggplot2, stringr, tidyr, patchwork, readr
###############################################################################

suppressPackageStartupMessages({
  library(h2o)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
  library(patchwork)
  library(readr)
})

# -------------------------- USER: EDIT THESE PATHS ----------------------------
species_dirs <- list(
  E3 = list(
    B = "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/H2o_Results/B/thinned_E3B",
    A = "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/H2o_Results/A/thinned_E3A"
  ),
  E4 = list(
    B = "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/H2o_Results/B/thinned_E4B",
    A = "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/H2o_Results/A/thinned_E4A"
  ),
  E5 = list(
    B = "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/H2o_Results/B/thinned_E5B",
    A = "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/H2o_Results/A/thinned_E5A"
  )
)

out_dir <- "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/VarImp_Simple_ALL_labels"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------ H2O -------------------------------------------
h2o.init()

# ------------------------------ LABELING --------------------------------------
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

# Length-stable, robust pretty-namer
pretty_var <- function(x) {
  x_chr  <- as.character(x)
  lower  <- tolower(x_chr)
  
  # Start with originals
  out <- x_chr
  
  # 1) BIOCLIM direct map (position-based; length-stable)
  map_idx <- match(lower, names(BIO_LUT))
  mapped  <- BIO_LUT[map_idx]
  repl    <- !is.na(map_idx)
  out[repl] <- mapped[repl]
  
  # 2) NDVI/EVI with year: NDVI_2024 -> "NDVI (2024)"
  idx_year <- grepl("^(ndvi|evi)_?\\d{4}$", lower, perl = TRUE)
  if (any(idx_year)) {
    yr <- str_extract(x_chr[idx_year], "\\d{4}$")
    base <- toupper(sub("_?\\d{4}$", "", x_chr[idx_year], perl = TRUE))
    out[idx_year] <- paste0(base, " (", yr, ")")
  }
  
  # 3) Plain NDVI/EVI
  out[grepl("^ndvi$", lower)] <- "NDVI"
  out[grepl("^evi$",  lower)] <- "EVI"
  
  # 4) Common topography/land layers (regex, case-insensitive)
  out[grepl("^elev|^elevation$", lower)] <- "Elevation (m)"
  out[grepl("^slope$", lower)]           <- "Slope (°)"
  out[grepl("^aspect$", lower)]          <- "Cos_Aspect (°)"
  out[grepl("^land\\s*cover$", lower)]   <- "LandCover"
  out[grepl("^forest\\s*cover$", lower)] <- "ForestCover"
  
  # 5) Remaining: make human-readable (underscores->spaces; Title Case)
  still_raw <- out == x_chr
  if (any(still_raw)) {
    tmp <- gsub("_", " ", out[still_raw])
    tmp <- str_trim(tmp)
    tmp <- str_squish(tmp)
    tmp <- str_to_title(tmp)
    out[still_raw] <- tmp
  }
  
  out
}

# ------------------------------ helpers ---------------------------------------
pick_model <- function(dir_path) {
  if (!dir.exists(dir_path)) stop("Dir not found: ", dir_path)
  m <- list.files(dir_path,
                  pattern = "StackedEnsemble_AllModels|StackedEnsemble_BestOfFamily",
                  full.names = TRUE)
  if (!length(m)) {
    m <- list.files(dir_path,
                    pattern = "XGBoost|GBM|DRF|DeepLearning|GLM",
                    full.names = TRUE)
  }
  if (!length(m)) return(NA_character_)
  m[1]
}

combined_varimp <- function(model) {
  is_ens <- grepl("StackedEnsemble", model@algorithm, ignore.case = TRUE)
  if (!is_ens) {
    vi <- try(h2o.varimp(model), silent = TRUE)
    if (inherits(vi, "try-error") || is.null(vi)) return(tibble(variable="", rel_imp=0)[0,])
    as.data.frame(vi)[, c("variable","relative_importance")] %>%
      group_by(variable) %>%
      summarise(rel_imp = sum(relative_importance), .groups="drop") %>%
      mutate(rel_imp = 100 * rel_imp / sum(rel_imp))
  } else {
    base_ids <- model@model$base_models
    base_ids <- base_ids[!grepl("Metalearner", base_ids, ignore.case = TRUE)]
    if (!length(base_ids)) return(tibble(variable="", rel_imp=0)[0,])
    parts <- lapply(base_ids, function(bid) {
      bm <- h2o.getModel(bid)
      vi <- try(h2o.varimp(bm), silent = TRUE)
      if (inherits(vi, "try-error") || is.null(vi)) return(NULL)
      as.data.frame(vi)[, c("variable","relative_importance")]
    })
    parts <- Filter(Negate(is.null), parts)
    if (!length(parts)) return(tibble(variable="", rel_imp=0)[0,])
    bind_rows(parts) %>%
      group_by(variable) %>%
      summarise(rel_imp = sum(relative_importance), .groups="drop") %>%
      mutate(rel_imp = 100 * rel_imp / sum(rel_imp))
  }
}

# Plot ALL variables by default
col_before <- "#7A7A7A"  # grey
col_after  <- "#1F77B4"  # blue
plot_dumbbell <- function(wide_df, sp, top_k = Inf) {
  dd <- wide_df %>%
    mutate(maxval = pmax(Before, After, na.rm = TRUE)) %>%
    arrange(desc(maxval))
  
  if (is.finite(top_k)) dd <- dd %>% slice_head(n = min(top_k, nrow(dd)))
  
  dd$label <- pretty_var(dd$variable)
  dd$label <- factor(dd$label, levels = rev(dd$label))
  
  ggplot(dd) +
    geom_segment(aes(x = label, xend = label, y = Before, yend = After),
                 linewidth = 0.9, color = "grey65") +
    geom_point(aes(label, Before), size = 3.2, shape = 21,
               fill = col_before, color = "black", stroke = 0.25) +
    geom_point(aes(label, After),  size = 3.2, shape = 21,
               fill = col_after,  color = "black", stroke = 0.25) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = paste0("Variable Importance — ", sp, " (Before vs After)"),
         x = NULL, y = "Relative importance (%)") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 10.5)
    )
}

figure_height <- function(n_vars, base = 2.0, per_var = 0.33, min_h = 6.0, max_h = 20.0) {
  h <- base + per_var * n_vars
  max(min_h, min(h, max_h))
}

# ------------------------------ main ------------------------------------------
all_wide <- list()
all_plots <- list()

for (sp in names(species_dirs)) {
  cat("\n=== ", sp, " ===\n", sep = "")
  dir_B <- species_dirs[[sp]]$B
  dir_A <- species_dirs[[sp]]$A
  
  path_B <- pick_model(dir_B)
  path_A <- pick_model(dir_A)
  if (is.na(path_B) || is.na(path_A)) {
    cat("  ! Missing model (B or A). Skipping ", sp, ".\n", sep = "")
    next
  }
  
  cat("  Loading B: ", basename(path_B), "\n", sep = "")
  cat("  Loading A: ", basename(path_A), "\n", sep = "")
  mdl_B <- h2o.loadModel(path_B)
  mdl_A <- h2o.loadModel(path_A)
  
  vi_B <- combined_varimp(mdl_B) %>% rename(Before = rel_imp)
  vi_A <- combined_varimp(mdl_A) %>% rename(After  = rel_imp)
  
  wide <- full_join(vi_B, vi_A, by = "variable") %>%
    mutate(across(c(Before, After), ~ replace_na(., 0)),
           delta = After - Before)
  
  # Add pretty labels to CSV (length-stable)
  wide <- wide %>%
    mutate(variable_label = pretty_var(variable)) %>%
    arrange(desc(pmax(Before, After)))
  
  all_wide[[sp]] <- wide
  write_csv(wide, file.path(out_dir, paste0("VarImp_BeforeAfter_", sp, ".csv")))
  
  p <- plot_dumbbell(wide, sp, top_k = Inf)
  h_px <- figure_height(nrow(wide))
  ggsave(file.path(out_dir, paste0("Fig_VarImp_Dumbbell_", sp, ".png")),
         p, width = 9.2, height = h_px, dpi = 320)
  ggsave(file.path(out_dir, paste0("Fig_VarImp_Dumbbell_", sp, ".pdf")),
         p, width = 9.2, height = h_px, device = cairo_pdf)
  all_plots[[sp]] <- p
}

if (length(all_plots)) {
  panel <- wrap_plots(all_plots[names(all_plots) %in% c("E3","E4","E5")], ncol = 1)
  ggsave(file.path(out_dir, "Fig_VarImp_Dumbbell_All.png"),
         panel, width = 9.6, height = 20, dpi = 320)
  ggsave(file.path(out_dir, "Fig_VarImp_Dumbbell_All.pdf"),
         panel, width = 9.6, height = 20, device = cairo_pdf)
}

if (length(all_wide)) {
  long_all <- bind_rows(lapply(names(all_wide), function(sp) {
    all_wide[[sp]] %>% mutate(species = sp, .before = 1)
  }))
  write_csv(long_all, file.path(out_dir, "VarImp_BeforeAfter_ALL.csv"))
}

cat("\nSaved outputs to: ", normalizePath(out_dir), "\n", sep = "")
try(h2o.shutdown(prompt = FALSE), silent = TRUE)
###############################################################################
