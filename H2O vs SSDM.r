################################################################################
# Habitat Suitability Comparison — SSDM vs H2O AutoML (Paper A)
# improves robustness, alignment, and ready figures.
################################################################################

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(readr)
  library(tools)
  library(cowplot)
  library(gridExtra)
  library(scales)
  library(rlang)
})

set.seed(123)

# ----------------------------- USER PATHS -------------------------------------
ssdm_dir <- "D:/Harin/Projects/SDM/BTEH/Results/SSDM_all/rasters"
h2o_dir  <- "D:/Harin/Projects/SDM/BTEH/H2o_Results/rasters"
out_dir  <- "D:/Harin/Projects/SDM/BTEH/Results/PaperA_Temporal"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dir.create(file.path(out_dir, "01_between_methods/rasters"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "01_between_methods/plots"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "02_temporal/rasters"),        recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "02_temporal/plots"),          recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "03_maps"),                    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "04_panels"),                  recursive = TRUE, showWarnings = FALSE)

# README for coauthors
cat(
  "Folders:\n",
  " 01_between_methods: pixel-wise SSDM vs H2O (diff rasters + metrics plots)\n",
  " 02_temporal: Δ (A−B) rasters, GSL rasters, and threshold metrics\n",
  " 03_maps: base suitability maps per dataset & method\n",
  " 04_panels: one figure per (E3|E4|E5) × (SSDM|H2O)\n",
  file = file.path(out_dir, "README.txt")
)

# -------------------------- DISCOVERY & PARSING -------------------------------
parse_meta <- function(fp, method_hint){
  nm <- basename(fp)
  # Expect tokens like: E1B, E2B, E3A, E4A, E5A
  dataset <- stringr::str_match(nm, "(E[1-5][AB])")[,2]
  tibble(file = fp, dataset = dataset, method = method_hint) |> tidyr::drop_na(dataset)
}

ssdm_files <- list.files(ssdm_dir, pattern = "\\.tif$", full.names = TRUE)
h2o_files  <- list.files(h2o_dir,  pattern = "\\.tif$", full.names = TRUE)

meta_ssdm <- dplyr::bind_rows(lapply(ssdm_files, parse_meta, method_hint = "SSDM"))
meta_h2o  <- dplyr::bind_rows(lapply(h2o_files,  parse_meta, method_hint = "H2O"))

# Pair SSDM/H2O by dataset
meta_all <- dplyr::full_join(meta_ssdm, meta_h2o, by = "dataset", suffix = c("_ssdm", "_h2o"))
datasets_all <- meta_all$dataset
stopifnot(length(datasets_all) > 0)

# --------------------------- UTILS & HELPERS ----------------------------------
align_to <- function(r1, r2, categorical = FALSE) {
  # Project -> resample -> mask -> crop to ensure exact alignment & NA mask match
  if (!compareGeom(r1, r2, stopOnError = FALSE)) {
    r2 <- project(r2, crs(r1), method = if (categorical) "near" else "bilinear")
    r2 <- resample(r2, r1,      method = if (categorical) "near" else "bilinear")
    r2 <- crop(r2, r1)
  }
  # apply common NA mask (intersection)
  m <- !is.na(r1) & !is.na(r2)
  r1 <- mask(r1, m, maskvalues = 0)
  r2 <- mask(r2, m, maskvalues = 0)
  list(r1 = r1, r2 = r2)
}

r_to_df_full <- function(r) {
  as.data.frame(r, xy = TRUE, na.rm = TRUE) |>
    stats::setNames(c("x","y","val"))
}

plot_raster_continuous <- function(r, title, out_png, center0 = FALSE) {
  df <- r_to_df_full(r)
  p <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(hjust = 0))
  if (center0) {
    lim <- max(abs(range(df$val, na.rm = TRUE)))
    p <- p + scale_fill_gradientn(colors = c("#2c7fb8","#ffffbf","#d7191c"),
                                  limits = c(-lim, lim), na.value = NA)
  } else {
    p <- p + scale_fill_viridis_c(na.value = NA)
  }
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
}

plot_raster_discrete <- function(r, title, out_png) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x","y","val")
  df$val <- factor(df$val, levels = c(-1,0,1), labels = c("Loss","Stable","Gain"))
  p <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    coord_equal() +
    scale_fill_manual(values = c("Loss"="#b2182b","Stable"="#f7f7f7","Gain"="#2166ac"), drop = FALSE) +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          plot.title = element_text(hjust = 0))
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
}

jaccard_binary <- function(b1, b2) {
  inter <- global(b1 & b2, "sum", na.rm = TRUE)[[1]]
  union <- global(b1 | b2, "sum", na.rm = TRUE)[[1]]
  ifelse(union == 0, NA_real_, inter / union)
}

quick_bar <- function(df, metric, ylab, outpng){
  ord <- df$dataset %>% factor(levels = sort(unique(df$dataset), decreasing = TRUE))
  p <- ggplot(df, aes(x = ord, y = .data[[metric]])) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = ylab, title = paste("Between-method:", ylab)) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0))
  ggsave(outpng, p, width = 7.5, height = 5, dpi = 300)
}

# ------------------ 01) BETWEEN-METHOD COMPARISONS (ALL 8) -------------------
bm_metrics <- list()

for (ds in datasets_all) {
  f_ssdm <- meta_all %>% filter(dataset == ds) %>% pull(file_ssdm)
  f_h2o  <- meta_all %>% filter(dataset == ds) %>% pull(file_h2o)
  if (length(f_ssdm) * length(f_h2o) == 0 || is.na(f_ssdm) || is.na(f_h2o)) {
    message("Skipping ", ds, " (missing pair)."); next
  }
  message("Between-method: ", ds)
  
  r_ssdm <- rast(f_ssdm)
  r_h2o  <- rast(f_h2o)
  al     <- align_to(r_ssdm, r_h2o, categorical = FALSE)
  r_ssdm <- al$r1; r_h2o <- al$r2
  
  # Base maps
  plot_raster_continuous(r_ssdm, paste0(ds, " — SSDM suitability"),
                         file.path(out_dir, "03_maps", paste0(ds, "_SSDM.png")))
  plot_raster_continuous(r_h2o,  paste0(ds, " — H2O suitability"),
                         file.path(out_dir, "03_maps", paste0(ds, "_H2O.png")))
  
  # Metrics
  v1 <- values(r_ssdm, mat = FALSE)
  v2 <- values(r_h2o,  mat = FALSE)
  idx <- !is.na(v1) & !is.na(v2)
  pearson  <- suppressWarnings(cor(v1[idx], v2[idx], method = "pearson"))
  spearman <- suppressWarnings(cor(v1[idx], v2[idx], method = "spearman"))
  rmse     <- sqrt(mean((v1[idx] - v2[idx])^2))
  mae      <- mean(abs(v1[idx] - v2[idx]))
  
  # Hotspot Jaccard (top 25%)
  q1 <- quantile(v1[idx], 0.75, na.rm = TRUE)
  q2 <- quantile(v2[idx], 0.75, na.rm = TRUE)
  h1 <- r_ssdm > q1
  h2 <- r_h2o  > q2
  jac_hot <- jaccard_binary(h1, h2)
  
  # Difference raster H2O − SSDM
  diff_h2o_minus_ssdm <- r_h2o - r_ssdm
  writeRaster(diff_h2o_minus_ssdm,
              file.path(out_dir, "01_between_methods/rasters",
                        paste0("diff_H2O_minus_SSDM_", ds, ".tif")),
              overwrite = TRUE, gdal = c("COMPRESS=LZW","BIGTIFF=IF_SAFER"), datatype = "FLT4S")
  plot_raster_continuous(diff_h2o_minus_ssdm,
                         paste0("Between-method difference (H2O − SSDM): ", ds),
                         file.path(out_dir, "01_between_methods/plots",
                                   paste0("diff_H2O_minus_SSDM_", ds, ".png")),
                         center0 = TRUE)
  
  bm_metrics[[ds]] <- tibble(
    dataset = ds, pearson = pearson, spearman = spearman,
    rmse = rmse, mae = mae, jaccard_hotspot_q75 = jac_hot
  )
}

bm_metrics_df <- bind_rows(bm_metrics)
write_csv(bm_metrics_df, file.path(out_dir, "01_between_methods", "between_method_metrics.csv"))

# Quick diagnostic bars
quick_bar(bm_metrics_df, "pearson", "Pearson r",
          file.path(out_dir, "01_between_methods/plots", "pearson_bar.png"))
quick_bar(bm_metrics_df, "rmse", "RMSE",
          file.path(out_dir, "01_between_methods/plots", "rmse_bar.png"))
quick_bar(bm_metrics_df, "jaccard_hotspot_q75", "Hotspot Jaccard (Q75)",
          file.path(out_dir, "01_between_methods/plots", "jaccard_q75_bar.png"))

# ------------------ 02) TEMPORAL CHANGE  Δ = A − B  (E3, E4, E5) -------------
temporal_ids <- c("E3", "E4", "E5")
thresholds   <- c(0.25, 0.5, 0.75)

get_file <- function(id, period, method){
  ds  <- paste0(id, period)  # "A"=After, "B"=Before
  col <- ifelse(method == "SSDM", "file_ssdm", "file_h2o")
  val <- meta_all %>% filter(dataset == ds) %>% pull(all_of(col))
  if (length(val) == 0) NA_character_ else val
}

area_rows <- list()
jacc_rows <- list()

for (id in temporal_ids) {
  for (method in c("SSDM", "H2O")) {
    fB <- get_file(id, "B", method)  # BEFORE
    fA <- get_file(id, "A", method)  # AFTER
    if (is.na(fB) || is.na(fA) || !file.exists(fB) || !file.exists(fA)) {
      message("Temporal Δ: skipping ", id, " | ", method, " (missing A or B)"); next
    }
    
    message("Temporal Δ (A − B): ", id, " | ", method)
    rB <- rast(fB); rA <- rast(fA)
    al <- align_to(rA, rB, categorical = FALSE)  # align B to A’s grid & mask
    rA <- al$r1; rB <- al$r2
    
    # Δ suitability = After − Before
    rDelta <- rA - rB
    writeRaster(rDelta,
                file.path(out_dir, "02_temporal/rasters",
                          paste0("delta_A_minus_B_", id, "_", method, ".tif")),
                overwrite = TRUE, gdal = c("COMPRESS=LZW","BIGTIFF=IF_SAFER"), datatype = "FLT4S")
    plot_raster_continuous(rDelta,
                           paste0("Temporal change Δ (A − B): ", id, " — ", method),
                           file.path(out_dir, "02_temporal/plots",
                                     paste0("delta_A_minus_B_", id, "_", method, ".png")),
                           center0 = TRUE)
    
    # Threshold-based Gain/Stable/Loss and metrics
    for (thr in thresholds) {
      bA <- rA >= thr   # AFTER
      bB <- rB >= thr   # BEFORE
      
      # Align categorical exactly (nearest) and common mask
      alB <- align_to(bA, bB, categorical = TRUE)
      bA  <- alB$r1; bB <- alB$r2
      
      # Jaccard(A,B)
      J <- jaccard_binary(bA, bB)
      
      # Areas (count of TRUE pixels)
      A_area <- global(bA, "sum", na.rm = TRUE)[[1]]
      B_area <- global(bB, "sum", na.rm = TRUE)[[1]]
      I_area <- global(bA & bB, "sum", na.rm = TRUE)[[1]]
      
      loss   <- (B_area - I_area)           # was suitable, now not
      gain   <- (A_area - I_area)           # newly suitable
      stable <- I_area
      pct_change <- ifelse(B_area == 0, NA_real_, (A_area - B_area) / B_area * 100)
      
      # GSL map: -1=Loss, 0=Stable, +1=Gain
      loss_r   <- (bB & (!bA))
      gain_r   <- ((!bB) & bA)
      stable_r <- (bA & bB)
      gsl <- (gain_r*1) + (stable_r*0) + (loss_r*(-1))
      names(gsl) <- "GSL"
      
      out_gsl_tif <- file.path(out_dir, "02_temporal/rasters",
                               paste0("GSL_", id, "_", method, "_thr", format(thr, nsmall = 2), ".tif"))
      writeRaster(gsl, out_gsl_tif, overwrite = TRUE,
                  gdal = c("COMPRESS=LZW","BIGTIFF=IF_SAFER"), datatype = "INT1S")
      plot_raster_discrete(gsl,
                           sprintf("Gain / Stable / Loss (thr=%.2f): %s — %s", thr, id, method),
                           file.path(out_dir, "02_temporal/plots",
                                     paste0("GSL_", id, "_", method, "_thr", thr, ".png"))
      )
      
      jacc_rows[[length(jacc_rows)+1]] <- tibble(id = id, method = method, threshold = thr, jaccard = J)
      area_rows[[length(area_rows)+1]] <- tibble(id = id, method = method, threshold = thr,
                                                 B_area = B_area, A_area = A_area,
                                                 loss = loss, stable = stable, gain = gain,
                                                 pct_change = pct_change)
    }
    
    # Histogram of Δ
    dfDelta <- r_to_df_full(rDelta)
    p_hist <- ggplot(dfDelta, aes(x = val)) +
      geom_histogram(bins = 60) +
      labs(title = paste0("Δ (A − B) distribution: ", id, " — ", method),
           x = "Suitability change", y = "Frequency") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(hjust = 0))
    ggsave(file.path(out_dir, "02_temporal/plots",
                     paste0("delta_hist_", id, "_", method, ".png")),
           p_hist, width = 7, height = 5, dpi = 300)
  }
}

jacc_df <- bind_rows(jacc_rows)
area_df <- bind_rows(area_rows)
write_csv(jacc_df, file.path(out_dir, "02_temporal", "jaccard_by_threshold.csv"))
write_csv(area_df, file.path(out_dir, "02_temporal", "area_change_by_threshold.csv"))

# Plots for temporal metrics
p_jacc <- ggplot(jacc_df, aes(x = factor(threshold), y = jaccard, group = method, linetype = method)) +
  geom_point() + geom_line() +
  facet_wrap(~ id, ncol = 3) +
  labs(x = "Threshold", y = "Jaccard(A,B)", title = "Temporal overlap across thresholds") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0))
ggsave(file.path(out_dir, "02_temporal/plots", "jaccard_lines.png"), p_jacc, width = 8, height = 4.8, dpi = 300)

p_area <- ggplot(area_df, aes(x = factor(threshold), y = pct_change, fill = method)) +
  geom_col(position = position_dodge(width = 0.75)) +
  facet_wrap(~ id, ncol = 3) +
  labs(x = "Threshold", y = "% area change (A vs B)", title = "Temporal area change by threshold") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(hjust = 0))
ggsave(file.path(out_dir, "02_temporal/plots", "area_change_bars.png"), p_area, width = 8.2, height = 4.8, dpi = 300)

# ==================================================================================
# 04) PANEL COMPOSER — one PNG per elephant × method
# ==================================================================================
panel_dir <- file.path(out_dir, "04_panels")

read_r <- function(fp) rast(fp)
exists_file <- function(x) is.character(x) && nzchar(x) && file.exists(x)

plot_gsl_mini <- function(r, title) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x","y","val")
  df$val <- factor(df$val, levels = c(-1,0,1), labels = c("Loss","Stable","Gain"))
  ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    coord_equal() +
    scale_fill_manual(values = c("Loss"="#b2182b","Stable"="#f7f7f7","Gain"="#2166ac"), drop = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          legend.position = "none", plot.title = element_text(hjust = 0))
}

plot_cont_mini <- function(r, title, center0 = FALSE) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x","y","val")
  p <- ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_raster() +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          legend.position = "none", plot.title = element_text(hjust = 0))
  if (center0) {
    lim <- max(abs(range(df$val, na.rm = TRUE)))
    p + scale_fill_gradientn(colors = c("#2c7fb8","#ffffbf","#d7191c"),
                             limits = c(-lim, lim), na.value = NA)
  } else {
    p + scale_fill_viridis_c(na.value = NA)
  }
}

# Load previously saved temporal CSVs
jacc_df <- read_csv(file.path(out_dir, "02_temporal", "jaccard_by_threshold.csv"), show_col_types = FALSE)
area_df <- read_csv(file.path(out_dir, "02_temporal", "area_change_by_threshold.csv"), show_col_types = FALSE)

# Build compact summary wide->long
summ_tbl <- jacc_df %>%
  mutate(threshold = as.character(threshold)) %>%
  pivot_wider(id_cols = c(id, method),
              names_from = threshold, values_from = jaccard,
              names_prefix = "Jacc_") %>%
  left_join(
    area_df %>%
      select(id, method, threshold, pct_change) %>%
      mutate(threshold = as.character(threshold)) %>%
      pivot_wider(id_cols = c(id, method),
                  names_from = threshold, values_from = pct_change,
                  names_prefix = "PctΔ_"),
    by = c("id", "method")
  ) %>%
  mutate(across(starts_with("Jacc_"), ~round(.x, 3)),
         across(starts_with("PctΔ_"), ~round(.x, 1)))

compose_panel <- function(id, method) {
  fA <- meta_all %>% filter(dataset == paste0(id, "A")) %>%
    pull(ifelse(method == "SSDM","file_ssdm","file_h2o")) %>% .[1]
  fB <- meta_all %>% filter(dataset == paste0(id, "B")) %>%
    pull(ifelse(method == "SSDM","file_ssdm","file_h2o")) %>% .[1]
  if (!exists_file(fA) || !exists_file(fB)) { message("Panel: missing A/B for ", id, " | ", method); return(invisible(NULL)) }
  
  rA <- read_r(fA); rB <- read_r(fB)
  al <- align_to(rA, rB, categorical = FALSE)  # align B to A’s grid & mask
  rA <- al$r1; rB <- al$r2
  
  fDelta <- file.path(out_dir, "02_temporal/rasters",
                      paste0("delta_A_minus_B_", id, "_", method, ".tif"))
  if (!exists_file(fDelta)) {
    writeRaster(rA - rB, fDelta, overwrite = TRUE,
                gdal = c("COMPRESS=LZW","BIGTIFF=IF_SAFER"), datatype = "FLT4S")
  }
  rDelta <- read_r(fDelta)
  
  thr_vals <- c(0.25, 0.5, 0.75)
  gsl_files <- sprintf(file.path(out_dir, "02_temporal/rasters",
                                 "GSL_%s_%s_thr%.2f.tif"), id, method, thr_vals)
  if (!all(file.exists(gsl_files))) { message("Panel: missing GSL for ", id, " | ", method); return(invisible(NULL)) }
  gsl_list <- lapply(seq_along(gsl_files), function(i) {
    r <- read_r(gsl_files[i])
    plot_gsl_mini(r, paste0("GSL (thr=", thr_vals[i], ")"))
  })
  
  pA     <- plot_cont_mini(rA,     paste0(id, " — ", method, " (After)"))
  pB     <- plot_cont_mini(rB,     paste0(id, " — ", method, " (Before)"))
  pDelta <- plot_cont_mini(rDelta, paste0("Δ (A − B)"), center0 = TRUE)
  
  mt <- summ_tbl %>% filter(id == !!id, method == !!method)
  if (nrow(mt) == 0) {
    mt_long <- tibble(
      Metric = c("Jaccard (thr=0.25)", "Jaccard (thr=0.50)", "Jaccard (thr=0.75)",
                 "%Δ area (thr=0.25)", "%Δ area (thr=0.50)", "%Δ area (thr=0.75)"),
      Value  = rep(NA_real_, 6)
    )
  } else {
    mt_long <- mt %>%
      select(starts_with("Jacc_"), starts_with("PctΔ_")) %>%
      pivot_longer(everything(), names_to = "MetricRaw", values_to = "Value") %>%
      mutate(Metric = case_when(
        MetricRaw == "Jacc_0.25" ~ "Jaccard (thr=0.25)",
        MetricRaw == "Jacc_0.5"  ~ "Jaccard (thr=0.50)",
        MetricRaw == "Jacc_0.75" ~ "Jaccard (thr=0.75)",
        MetricRaw == "PctΔ_0.25" ~ "%Δ area (thr=0.25)",
        MetricRaw == "PctΔ_0.5"  ~ "%Δ area (thr=0.50)",
        MetricRaw == "PctΔ_0.75" ~ "%Δ area (thr=0.75)",
        TRUE ~ MetricRaw
      )) %>% select(Metric, Value)
  }
  
  mt_long <- mt_long %>%
    mutate(Value_str = ifelse(is.na(Value), "NA",
                              ifelse(grepl("^Jaccard", Metric),
                                     sprintf("%.3f", Value),
                                     sprintf("%.1f", Value))))
  
  jacc_pal <- scales::col_numeric(palette = c("#e9f7e9", "#40a040", "#005a00"),
                                  domain = c(0, 1), na.color = "#f0f0f0")
  pct_vals <- mt_long %>% filter(grepl("^%Δ area", Metric)) %>% pull(Value)
  rng <- range(pct_vals, na.rm = TRUE)
  max_abs <- ifelse(all(is.finite(rng)), max(abs(rng)), 50)
  pct_pal <- scales::col_numeric(palette = c("#2166ac", "#ffffff", "#b2182b"),
                                 domain = c(-max_abs, max_abs), na.color = "#f0f0f0")
  
  row_bg <- vapply(seq_len(nrow(mt_long)), function(i) {
    if (grepl("^Jaccard", mt_long$Metric[i])) jacc_pal(mt_long$Value[i]) else pct_pal(mt_long$Value[i])
  }, FUN.VALUE = character(1))
  
  tbl_theme <- gridExtra::ttheme_minimal(
    core = list(
      fg_params = list(cex = 1.05, fontface = "plain", col = "black"),
      bg_params = list(fill = row_bg, col = NA)
    ),
    colhead = list(fg_params = list(cex = 1.15, fontface = "bold"))
  )
  
  tbl_grob <- gridExtra::tableGrob(
    mt_long %>% transmute(Metric, Value = Value_str),
    rows = NULL, theme = tbl_theme
  )
  tbl_grob$widths  <- grid::unit(c(6, 3), "cm")
  tbl_grob$heights <- grid::unit(rep(0.7, nrow(mt_long) + 1), "cm")
  
  top_row    <- plot_grid(pB, pA, pDelta, nrow = 1, rel_widths = c(1,1,1.1))
  gsl_strip  <- plot_grid(plotlist = gsl_list, nrow = 1)
  bottom_row <- plot_grid(gsl_strip, ggplot() + theme_void(), ggplot() + theme_void(),
                          tbl_grob, nrow = 1, rel_widths = c(2.3, 0.05, 0.05, 1.1))
  
  title <- ggdraw() +
    draw_label(paste0("Elephant ", id, " — ", method,
                      " | Before vs After: base maps, Δ (A−B), GSL strip, metrics"),
               x = 0, hjust = 0, fontface = "bold", size = 12) +
    theme(plot.margin = margin(5,5,5,5))
  
  panel <- plot_grid(title, top_row, bottom_row, ncol = 1, rel_heights = c(0.12, 1, 0.75))
  out_png <- file.path(panel_dir, paste0(id, "_", method, "_PANEL.png"))
  ggsave(out_png, panel, width = 14, height = 9, dpi = 300)
  message("Wrote: ", out_png)
}

for (id in c("E3","E4","E5")) {
  for (method in c("SSDM","H2O")) {
    compose_panel(id, method)
  }
}

message("All done. Outputs saved under: ", out_dir)
################################################################################
