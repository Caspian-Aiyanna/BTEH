suppressPackageStartupMessages({
  if (!requireNamespace("h2o", quietly = TRUE)) {
    stop("The 'h2o' package is required. Please install it before running the H2O scripts.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("The 'dplyr' package is required. Please install it before running the H2O scripts.")
  }
})

restrict_to_model_family <- function(leaderboard, preferred = c("StackedEnsemble", "GBM", "XGBoost", "DRF", "DeepLearning", "GLM")) {
  leaderboard <- dplyr::mutate(leaderboard, model_category = dplyr::case_when(
    grepl("StackedEnsemble", model_id) ~ "StackedEnsemble",
    grepl("GBM", model_id) ~ "GBM",
    grepl("XGBoost", model_id) ~ "XGBoost",
    grepl("DeepLearning", model_id) ~ "DeepLearning",
    grepl("DRF", model_id) ~ "DRF",
    grepl("GLM", model_id) ~ "GLM",
    TRUE ~ "Other"
  ))
  leaderboard <- leaderboard[order(match(leaderboard$model_category, preferred), leaderboard$rank), ]
  leaderboard[!is.na(leaderboard$model_category) & leaderboard$model_category != "Other", , drop = FALSE]
}

predict_raster_h2o <- function(env_stack, leader, out_path, batch = 200000) {
  bs <- terra::blockSize(env_stack, minblocks = max(1, terra::nlyr(env_stack)))
  out <- terra::rast(env_stack[[1]])
  names(out) <- "suitability"
  terra::writeStart(out, out_path, overwrite = TRUE)
  on.exit(try(terra::writeStop(out), silent = TRUE), add = TRUE)
  for (i in seq_len(bs$n)) {
    vals <- terra::getValues(env_stack, row = bs$row[i], nrows = bs$nrows[i], mat = TRUE)
    n <- nrow(vals)
    if (is.null(n) || n == 0) {
      terra::writeValues(out, numeric(0), bs$row[i])
      next
    }
    pred_block <- rep(NA_real_, n)
    idx_all <- which(rowSums(is.na(vals)) == 0)
    if (length(idx_all) > 0) {
      for (start in seq(1, length(idx_all), by = batch)) {
        idx <- idx_all[start:min(start + batch - 1, length(idx_all))]
        dfb <- as.data.frame(vals[idx, , drop = FALSE])
        hf <- h2o::as.h2o(dfb)
        p <- h2o::h2o.predict(leader, hf)[["p1"]]
        pred_block[idx] <- as.vector(p)
      }
    }
    terra::writeValues(out, pred_block, bs$row[i])
  }
  terra::writeStop(out)
  invisible(out_path)
}
