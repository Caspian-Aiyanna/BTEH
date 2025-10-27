# BTEH Species Distribution Modeling Pipeline

This repository provides a modular, reproducible workflow for species distribution
modeling (SDM) using H2O AutoML and SSDM comparisons. The project layout is fully
described in `config.yml` and is driven via a simple `Makefile` so that each stage
can be resumed or re-run independently without modifying the underlying modeling
logic.

## Project Layout

```
BTEH/
├─ renv.lock
├─ .Rprofile
├─ config.yml
├─ data/
│  ├─ raw/
│  ├─ clean/
│  ├─ occ/
│  └─ envi/
├─ plans/
│  ├─ A/
│  └─ B/
├─ results/
│  ├─ H2O/
│  ├─ SSDM/
│  └─ compare/
├─ logs/
├─ R/
│  ├─ utils_io.R
│  ├─ utils_repro.R
│  ├─ utils_kendall.R
│  ├─ utils_h2o.R
│  └─ utils_plot.R
├─ scripts/
│  ├─ 02_dbscan_thin_degrees.R
│  ├─ 03_h2o_train.R
│  ├─ 05_h20_vs_ssdm_results.R
│  ├─ 05_uncertainity.R
│  └─ variable_processing_app.R
└─ Makefile
```

Each dataset receives its own sub-folders in `data/`, `plans/`, and `results/`.
Deterministic artifacts such as the Kendall correlation keep/drop lists are
stored in `plans/<dataset>/` so that the workflow can be resumed safely.

## Configuration

All tunable parameters live in `config.yml`. Key sections include:

- `project`: global seed and reproducibility settings.
- `paths`: relative locations for data, plans, results, and logs.
- `thinning`: DBSCAN parameters used by `02_dbscan_thin_degrees.R`.
- `modeling`: H2O AutoML budgets, spatial CV settings, and raster prediction
  batch sizes.
- `datasets`: dataset-specific occurrence files, environmental raster folders,
  and output destinations.
- `uncertainty`: before/after model directories for variable importance
  comparisons.

Update `config.yml` when adding new datasets or adjusting modeling knobs. The
utility functions in `R/` resolve all paths relative to the project root so the
pipeline remains portable.

## Running the Pipeline

The `Makefile` exposes the main steps:

```bash
# Thinning occurrence records for the configured dataset (defaults to runtime.dataset)
make dbscan

# Train H2O AutoML models for the dataset defined in config.yml
make h2o

# Generate SSDM vs H2O comparison figures/tables
make compare

# Build before/after variable-importance panels
make uncertainty

# Execute the full chain
make all
```

Set `DATASET` (e.g. `make DATASET=B h2o`) to override the dataset at runtime.

## Shiny Raster Processing App

`scripts/variable_processing_app.R` provides a Shiny dashboard for preparing
raster stacks with uniform resolution, extent, and CRS. Launch it with
`Rscript scripts/variable_processing_app.R` or run interactively from RStudio.

## Reproducibility

`.Rprofile` automatically activates the `renv` project if installed. Use
`renv::restore()` to install the frozen dependencies captured in `renv.lock`.

## Logging & Outputs

Intermediate artifacts (e.g., Kendall keep/drop lists, AutoML leaderboards,
comparison metrics) are saved under the dataset-specific folders described in
`config.yml`. Logs can be placed in the `logs/` directory as desired.

