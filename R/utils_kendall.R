# Kendall correlation utilities shared across scripts

prune_by_kendall <- function(cmat, cutoff = 0.8) {
  M <- abs(cmat)
  diag(M) <- 0
  if (ncol(M) < 2) {
    return(character(0))
  }
  to_drop <- character(0)
  keep <- colnames(M)
  repeat {
    maxv <- suppressWarnings(max(M, na.rm = TRUE))
    if (!is.finite(maxv) || maxv < cutoff) {
      break
    }
    idx <- which(M == maxv, arr.ind = TRUE)[1, ]
    c1 <- colnames(M)[idx[1]]
    c2 <- colnames(M)[idx[2]]
    m1 <- mean(M[idx[1], -idx[1]], na.rm = TRUE)
    m2 <- mean(M[idx[2], -idx[2]], na.rm = TRUE)
    drop_var <- if (is.na(m1) || is.na(m2)) c2 else if (m1 >= m2) c1 else c2
    to_drop <- c(to_drop, drop_var)
    keep <- setdiff(keep, drop_var)
    if (length(keep) < 2) {
      break
    }
    M <- M[keep, keep, drop = FALSE]
  }
  unique(to_drop)
}

compute_kendall_matrix <- function(env_stack, sample_size = 50000) {
  nsamp <- min(sample_size, terra::ncell(env_stack[[1]]))
  samp_pts <- terra::spatSample(env_stack[[1]], size = nsamp, method = "random",
                                as.points = TRUE, na.rm = TRUE)
  samp_vals <- terra::extract(env_stack, samp_pts)[, -1, drop = FALSE]
  samp_df <- as.data.frame(samp_vals)
  samp_df <- stats::na.omit(samp_df)
  if (nrow(samp_df) > sample_size) {
    samp_df <- samp_df[sample.int(nrow(samp_df), sample_size), , drop = FALSE]
  }
  suppressWarnings(stats::cor(samp_df, method = "kendall",
                              use = "pairwise.complete.obs"))
}

