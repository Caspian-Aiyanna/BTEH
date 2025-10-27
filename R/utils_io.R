# Utility functions for loading configuration and managing project paths

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The 'yaml' package is required. Please install it before running the pipeline.")
  }
  if (!requireNamespace("fs", quietly = TRUE)) {
    stop("The 'fs' package is required. Please install it before running the pipeline.")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("The 'readr' package is required. Please install it before running the pipeline.")
  }
})

#' Load project configuration
#'
#' @param config_path Path to the YAML configuration file relative to the
#'   project root (defaults to `config.yml`).
#' @return A list representing the parsed configuration with the project root
#'   stored in the `__root` element.
get_config <- function(config_path = "config.yml") {
  cfg <- yaml::read_yaml(config_path)
  cfg$`__root` <- fs::path_abs(".")
  class(cfg) <- c("bteh_config", class(cfg))
  cfg
}

#' Resolve a path relative to the project root
#'
#' @param cfg Configuration object returned by `get_config()`.
#' @param ... Path components to be joined.
#' @return Absolute path constructed from the project root and the provided
#'   components.
resolve_path <- function(cfg, ...) {
  stopifnot(inherits(cfg, "bteh_config"))
  fs::path_abs(fs::path(cfg$`__root`, ...))
}

#' Create a directory if it does not yet exist
#'
#' @param path Path to create.
#' @param recurse Whether to create intermediate directories.
#' @return The input path (invisibly).
ensure_dir <- function(path, recurse = TRUE) {
  if (!fs::dir_exists(path)) {
    fs::dir_create(path, recurse = recurse)
  }
  invisible(path)
}

#' Load a CSV artifact from the deterministic plan directory
#'
#' @param cfg Project configuration.
#' @param dataset Dataset identifier (e.g. "A").
#' @param filename Name of the CSV file inside the dataset plan directory.
#' @param required Logical; if `TRUE` (default) the file must exist.
#' @return Tibble with the file contents or `NULL` if not required and missing.
read_plan_csv <- function(cfg, dataset, filename, required = TRUE) {
  plan_dir <- resolve_path(cfg, cfg$paths$plans, dataset)
  file <- fs::path(plan_dir, filename)
  if (!fs::file_exists(file)) {
    if (required) {
      stop(sprintf("Missing plan artifact: %s", file))
    }
    return(NULL)
  }
  readr::read_csv(file, show_col_types = FALSE)
}

#' Write a CSV artifact into the deterministic plan directory
#'
#' @param cfg Project configuration.
#' @param dataset Dataset identifier (e.g. "A").
#' @param filename File name relative to the dataset plan directory.
#' @param data Data frame to write.
#' @return Invisibly, the written file path.
write_plan_csv <- function(cfg, dataset, filename, data) {
  plan_dir <- resolve_path(cfg, cfg$paths$plans, dataset)
  ensure_dir(plan_dir)
  file <- fs::path(plan_dir, filename)
  readr::write_csv(data, file)
  invisible(file)
}

#' Safely read a CSV file (wrapper around readr::read_csv)
#'
#' @param path Path to the CSV file.
#' @param ... Additional arguments passed to `readr::read_csv()`.
#' @return Tibble with the file contents.
read_csv_safely <- function(path, ...) {
  readr::read_csv(path, show_col_types = FALSE, ...)
}

#' Safely write a CSV file creating parent directories when needed
#'
#' @param data Data frame to write.
#' @param path Output path.
#' @param ... Additional arguments passed to `readr::write_csv()`.
write_csv_safely <- function(data, path, ...) {
  ensure_dir(fs::path_dir(path))
  readr::write_csv(data, path, ...)
  invisible(path)
}

