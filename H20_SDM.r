################################################################################
# Paper A — Single Dataset Runner (Current-Day, 2.5 km spatial CV)
# - Kendall’s tau collinearity screen (|τ| >= 0.8)
# - H2O AutoML on full data (leader saved)
# - 5-fold spatial CV with 2.5 km blocks (AUC/RMSE/logloss)
# - Variable importance + partial plots
# - Memory-safe chunked GeoTIFF prediction
#
# Requires: terra, h2o, dplyr, readr, ggplot2, tidyr, sf, blockCV, tools
################################################################################

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
})

set.seed(123)

# ----------------------------- USER PATHS -------------------------------------
# Set these 3 paths for ONE dataset
env_dir <- "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/envi/B"      # folder of .tif env variables
csv_path <- "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Balu/Balu_rep1.csv"  # presence CSV with lon,lat
res_dir <- "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Balu/output"              # output folder

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

# Tunables
block_km        <- 2.5
n_folds         <- 5
kendall_cutoff  <- 0.8
kendall_sample  <- 50000L
partial_vars_n  <- 5
automl_models   <- 30       # full-data AutoML budget
cv_automl_models <- 25      # per-fold budget
prediction_batch <- 200000  # number of cells per h2o batch during raster prediction

# --------------------- helpers: Kendall prune + chunked predict ---------------
prune_by_kendall <- function(cmat, cutoff = 0.8) {
  M <- abs(cmat); diag(M) <- 0
  if (ncol(M) < 2) return(character(0))
  to_drop <- character(0); keep <- colnames(M)
  repeat {
    maxv <- suppressWarnings(max(M, na.rm = TRUE))
    if (!is.finite(maxv) || maxv < cutoff) break
    idx <- which(M == maxv, arr.ind = TRUE)[1, ]
    c1  <- colnames(M)[idx[1]]; c2 <- colnames(M)[idx[2]]
    m1  <- mean(M[idx[1], -idx[1]], na.rm = TRUE)
    m2  <- mean(M[idx[2], -idx[2]], na.rm = TRUE)
    drop_var <- if (is.na(m1) || is.na(m2)) c2 else if (m1 >= m2) c1 else c2
    to_drop  <- c(to_drop, drop_var)
    keep     <- setdiff(keep, drop_var)
    if (length(keep) < 2) break
    M <- M[keep, keep, drop = FALSE]
  }
  unique(to_drop)
}

# Chunked raster prediction with H2O model
predict_raster_h2o <- function(env_stack, leader, out_path, batch = 200000) {
  bs <- terra::blockSize(env_stack, minblocks = max(1, nlyr(env_stack)))
  out <- terra::rast(env_stack[[1]]); names(out) <- "suitability"
  w <- terra::writeStart(out, out_path, overwrite = TRUE)
  on.exit(try(terra::writeStop(out), silent = TRUE), add = TRUE)
  
  for (i in seq_len(bs$n)) {
    vals <- terra::getValues(env_stack, row = bs$row[i], nrows = bs$nrows[i], mat = TRUE)
    # vals is matrix [nrow_block x nvars]; keep NAs
    n <- nrow(vals); if (is.null(n) || n == 0) { terra::writeValues(out, rep(NA_real_, 0), bs$row[i]); next }
    # Split into batches to limit H2O memory
    pred_block <- rep(NA_real_, n)
    idx_all <- which(rowSums(is.na(vals)) == 0)
    if (length(idx_all) > 0) {
      for (start in seq(1, length(idx_all), by = batch)) {
        idx <- idx_all[start:min(start + batch - 1, length(idx_all))]
        dfb <- as.data.frame(vals[idx, , drop = FALSE])
        hf  <- as.h2o(dfb)
        p   <- h2o.predict(leader, hf)[["p1"]]
        pred_block[idx] <- as.vector(p)
      }
    }
    terra::writeValues(out, pred_block, bs$row[i])
    cat(sprintf("Predicted block %d/%d\n", i, bs$n))
  }
  terra::writeStop(out)
  invisible(out_path)
}

# --------------------------- H2O SESSION --------------------------------------
h2o.init()  # add max_mem_size="10G" if needed

# ---------------------- Load environmental rasters ----------------------------
env_files <- list.files(env_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(env_files) > 0)
env_names <- make.names(file_path_sans_ext(basename(env_files)), unique = TRUE)
env_stack <- rast(env_files); names(env_stack) <- env_names
all_vars  <- names(env_stack)

message("====================================================")
message("Single dataset run")
message("Rasters: ", length(env_files), " | Variables: ", paste0(head(all_vars, 5), ifelse(length(all_vars)>5, " ...", "")))

# ------------------ Kendall’s tau screening for collinearity ------------------
nsamp <- min(kendall_sample, ncell(env_stack[[1]]))
samp_pts  <- spatSample(env_stack[[1]], size = nsamp, method = "random", as.points = TRUE, na.rm = TRUE)
samp_vals <- terra::extract(env_stack, samp_pts)[, -1, drop = FALSE]
samp_df   <- as.data.frame(samp_vals)
samp_df   <- samp_df[stats::complete.cases(samp_df), , drop = FALSE]
if (nrow(samp_df) > kendall_sample) {
  samp_df <- samp_df[sample(nrow(samp_df), kendall_sample), , drop = FALSE]
}
message("Computing Kendall τ on ", nrow(samp_df), " rows × ", ncol(samp_df), " vars ...")
kt <- suppressWarnings(stats::cor(samp_df, method = "kendall", use = "pairwise.complete.obs"))

drop_vars <- prune_by_kendall(kt, cutoff = kendall_cutoff)
keep_vars <- setdiff(colnames(kt), drop_vars)

if (length(drop_vars) > 0) {
  message("Kendall filter (|τ| >= ", kendall_cutoff, "): dropping ", length(drop_vars), " vars: ",
          paste(utils::head(drop_vars, 10), collapse = ", "),
          ifelse(length(drop_vars) > 10, " ...", ""))
  env_stack <- env_stack[[keep_vars]]
  all_vars  <- keep_vars
} else {
  message("Kendall filter: no variables dropped at |τ| >= ", kendall_cutoff)
}

# Save heatmap + dropped list
try({
  kt_keep <- kt[keep_vars, keep_vars, drop = FALSE]
  kt_long <- as.data.frame(as.table(kt_keep))
  names(kt_long) <- c("Var1", "Var2", "tau")
  g <- ggplot(kt_long, aes(Var1, Var2, fill = tau)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, fill = "Kendall τ",
         title = paste0("Kendall correlation (kept variables, |τ|<", kendall_cutoff, ")")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(file.path(res_dir, "kendall_heatmap_kept.png"), g, width = 8, height = 6, dpi = 300)
  write.csv(data.frame(dropped_variable = drop_vars),
            file.path(res_dir, "kendall_dropped_vars.csv"), row.names = FALSE)
}, silent = TRUE)

# -------------------------- Load presence CSV ---------------------------------
sp <- read_csv(csv_path, show_col_types = FALSE)
stopifnot(all(c("lon","lat") %in% names(sp)))
dset <- file_path_sans_ext(basename(csv_path))
species_dir <- file.path(res_dir, dset)
dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
message("---- Training dataset: ", dset)

# Presence(1) + background(0) balanced
pres <- sp %>% dplyr::select(lon, lat) %>% mutate(pa = 1)
bg_pts    <- spatSample(env_stack[[1]], size = nrow(pres), method = "random", as.points = TRUE)
bg_coords <- crds(bg_pts)
bg_df     <- data.frame(lon = bg_coords[,1], lat = bg_coords[,2], pa = 0)

df_sp <- bind_rows(pres, bg_df)

# Extract raster values
pts  <- vect(df_sp, geom = c("lon","lat"), crs = crs(env_stack))
vals <- terra::extract(env_stack, pts)[,-1, drop = FALSE]
df   <- bind_cols(df_sp, as.data.frame(vals)) %>% tidyr::drop_na()

# ---------------- SPATIAL BLOCK CROSS-VALIDATION (2.5 km) --------------------
sf_pts <- st_as_sf(df, coords = c("lon","lat"), crs = crs(env_stack))

folds_obj <- NULL
if ("cv_spatial" %in% getNamespaceExports("blockCV")) {
  folds_obj <- blockCV::cv_spatial(
    x = sf_pts, column = "pa", size = block_km * 1000, k = n_folds
  )
  if (!is.null(folds_obj$folds_ids)) {
    df$fold_id <- as.integer(folds_obj$folds_ids)
  } else if (!is.null(folds_obj$folds_list)) {
    fold_id <- rep(NA_integer_, nrow(df))
    for (k in seq_along(folds_obj$folds_list)) fold_id[folds_obj$folds_list[[k]]] <- k
    df$fold_id <- fold_id
  } else stop("cv_spatial(): cannot find fold ids; check blockCV version.")
} else {
  sb <- blockCV::spatialBlock(
    speciesData = sf_pts, species = "pa", theRange = block_km * 1000,
    k = n_folds, selection = "systematic", biomod2Format = FALSE, showBlocks = FALSE
  )
  fold_id <- rep(NA_integer_, nrow(df))
  for (k in seq_along(sb$folds)) fold_id[sb$folds[[k]]$test] <- k
  df$fold_id <- fold_id
}
stopifnot(!any(is.na(df$fold_id)))

# ---------------- AutoML (FULL DATA) — Leader model saved ---------------------
hf <- as.h2o(df)
hf[,"pa"] <- h2o.asfactor(hf[,"pa"]) 
aml <- h2o.automl(x = all_vars, y = "pa", training_frame = hf,
                  max_models = automl_models, seed = 123)
leader <- aml@leader
model_path <- h2o.saveModel(leader, path = species_dir, force = TRUE)

# In-sample performance (reference; spatial CV below)
perf <- h2o.performance(leader, newdata = hf)
metrics <- data.frame(
  dataset = dset,
  rmse    = h2o.rmse(perf),
  auc     = h2o.auc(perf),
  logloss = h2o.logloss(perf)
)
write_csv(metrics, file.path(species_dir, paste0("metrics_in_sample_", dset, ".csv")))

# VarImp + partials
varimp <- tryCatch(as.data.frame(h2o.varimp(leader)), error = function(e) NULL)
if (!is.null(varimp) && nrow(varimp) > 0) {
  varimp$dataset <- dset
  write_csv(varimp, file.path(species_dir, paste0("varimp_", dset, ".csv")))
  pdf(file.path(species_dir, paste0("varimp_", dset, ".pdf")));  print(h2o.varimp_plot(leader, num_of_features = 10)); dev.off()
  pp_vars <- head(all_vars, partial_vars_n)
  suppressWarnings({
    pdf(file.path(species_dir, paste0("partial_", dset, ".pdf"))); print(h2o.partialPlot(object = leader, data = hf, cols = pp_vars)); dev.off()
  })
}

# ----------------------- Spatial CV over folds --------------------------------
cv_rows <- list()
for (k in sort(unique(df$fold_id))) {
  trn <- df[df$fold_id != k, , drop = FALSE]
  tst <- df[df$fold_id ==  k, , drop = FALSE]
  if (nrow(trn) < 50 || nrow(tst) < 20) next
  
  hf_trn <- as.h2o(trn);  hf_trn[,"pa"] <- h2o::asfactor(hf_trn[,"pa"])
  hf_tst <- as.h2o(tst);  hf_tst[,"pa"] <- h2o::asfactor(hf_tst[,"pa"])
  
  
  aml_cv <- h2o.automl(x = all_vars, y = "pa", training_frame = hf_trn,
                       max_models = cv_automl_models, seed = 123 + k)
  leader_cv <- aml_cv@leader
  perf_cv   <- h2o.performance(leader_cv, newdata = hf_tst)
  
  cv_rows[[length(cv_rows) + 1]] <- data.frame(
    dataset = dset, fold = k,
    auc     = h2o.auc(perf_cv),
    rmse    = h2o.rmse(perf_cv),
    logloss = h2o.logloss(perf_cv)
  )
}

if (length(cv_rows) > 0) {
  cv_df  <- bind_rows(cv_rows)
  cv_sum <- cv_df %>%
    summarise(
      dataset      = dset,
      folds        = n(),
      auc_mean     = mean(auc, na.rm = TRUE),
      auc_sd       = sd(auc, na.rm = TRUE),
      rmse_mean    = mean(rmse, na.rm = TRUE),
      rmse_sd      = sd(rmse, na.rm = TRUE),
      logloss_mean = mean(logloss, na.rm = TRUE),
      logloss_sd   = sd(logloss, na.rm = TRUE)
    )
  write_csv(cv_df,  file.path(species_dir, paste0("spatialCV_folds_", dset, ".csv")))
  write_csv(cv_sum, file.path(species_dir, paste0("spatialCV_summary_", dset, ".csv")))
  
  g1 <- ggplot(cv_df, aes(x = factor(fold), y = auc)) +
    geom_col() +
    labs(x = "Fold", y = "AUC (spatial CV)", title = paste("Spatial CV —", dset)) +
    theme_minimal(base_size = 11)
  ggsave(file.path(species_dir, paste0("spatialCV_auc_", dset, ".png")) , g1, width = 5.5, height = 3.6, dpi = 300)
}

# ---------------------- Save raster prediction (current) ----------------------
pred_tif <- file.path(species_dir, paste0("prediction_", dset, ".tif"))
predict_raster_h2o(env_stack, leader, pred_tif, batch = prediction_batch)
message("Saved prediction: ", pred_tif)

# ------------------------------ FINISH ----------------------------------------
try(h2o.shutdown(prompt = FALSE), silent = TRUE)
message("Done.")
################################################################################

