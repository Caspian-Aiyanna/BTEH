# Reproducibility helpers (seeds, parallel safety, logging)

safely_register_parallel <- function(single_core = TRUE) {
  if (!requireNamespace("future", quietly = TRUE)) {
    warning("future package not installed; parallel control skipped")
    return(invisible(FALSE))
  }
  if (single_core) {
    future::plan(future::sequential)
  } else {
    future::plan(future::multisession)
  }
  invisible(TRUE)
}

set_reproducible_seed <- function(seed) {
  if (!is.null(seed)) {
    set.seed(seed)
    if (requireNamespace("withr", quietly = TRUE)) {
      withr::local_seed(seed)
    }
  }
  invisible(seed)
}

with_repro_context <- function(expr, seed = NULL, single_core = TRUE) {
  set_reproducible_seed(seed)
  safely_register_parallel(single_core)
  force(expr)
}

