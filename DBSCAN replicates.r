###############################################################################
# 1. Setup & Data Loading
###############################################################################
# install.packages(c("dbscan", "dplyr", "sf"))

library(dbscan)
library(dplyr)
library(sf)

# Input directory (original datasets) and output directory
in_dir  <- "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Uncertainity/data/OG_data/thin"
out_dir <- "D:/PhD/Projects/SDM_projects/BTEH/SDM_ele/Balu"
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# List of the datasets
datasets   <- list.files(in_dir, pattern = ".csv$", full.names = TRUE)
data_list  <- lapply(datasets, read.csv)
base_names <- gsub(".csv$", "", basename(datasets))

###############################################################################
# 2. DBSCAN-based Thinning Function
###############################################################################
sample_with_dbscan <- function(df, 
                               eps = 0.004,     
                               minPts = 15,    
                               fraction = 0.09, 
                               min_samples = 4) {
  # Run DBSCAN on coordinates
  db <- dbscan(as.matrix(df[, c("lon", "lat")]), eps = eps, minPts = minPts)
  df$cluster_id <- db$cluster
  
  sampled_data <- df %>%
    group_by(cluster_id) %>%
    group_modify(~ {
      if (.y$cluster_id == 0) { 
        slice_sample(.x, prop = fraction)   # noise
      } else {  
        n_to_sample <- max(min_samples, round(nrow(.x) * fraction))
        slice_sample(.x, n = min(nrow(.x), n_to_sample))
      }
    }) %>%
    ungroup()
  
  return(sampled_data)
}

###############################################################################
# 3. Generate Three Replicates for Each Dataset
###############################################################################
eps_value       <- 0.004
minPts_value    <- 15
fraction_value  <- 0.09
min_samples_val <- 4

for (i in seq_along(data_list)) {
  df <- data_list[[i]]
  dataset_name <- base_names[i]
  
  for (rep in 1:3) {
    set.seed(rep)  # ensures different replicates but reproducible
    thinned <- sample_with_dbscan(
      df,
      eps = eps_value,
      minPts = minPts_value,
      fraction = fraction_value,
      min_samples = min_samples_val
    )
    
    # Output file name with replicate index
    out_file <- file.path(out_dir, paste0(dataset_name, "_rep", rep, ".csv"))
    write.csv(thinned, out_file, row.names = FALSE)
    
    cat("Saved:", out_file, "(", nrow(thinned), "points)\n")
  }
}
