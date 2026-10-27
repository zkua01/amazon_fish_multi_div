# =============================================================================
# Script:  02a_build_community_matrix.R
# Purpose: Build the site × species abundance matrix from cleaned occurrences
#          and imputed trait data.
#
# Input:   data/processed/amazon_occurrences_clean.csv
#          data/processed/trait_imputed_2492.csv
# Output:  output/intermediate/community_matrix.rds
# =============================================================================

library(tidyverse)
library(here)

source(here("scripts", "utils", "clean_name.R"))

input_dir  <- here("data", "processed")
output_dir <- here("output", "intermediate")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
abd        <- read.csv(file.path(input_dir, "amazon_occurrences_clean.csv"))
trait_fill <- read.csv(file.path(input_dir, "trait_imputed_2492.csv"))

# Trait matrix (species × traits)
trait <- trait_fill |>
  dplyr::select(Genus.species, MBl:CPt) |>
  mutate(Genus.species = clean_name(Genus.species)) |>
  tibble::column_to_rownames(var = "Genus.species")

# Community matrix (sites × species)
comm <- abd |>
  mutate(
    valid_name = clean_name(str_replace_all(valid_name, " ", "_")),
    comm_code  = paste0("Comm_", formatC(SubDrainage.Code, width = 3, format = "d", flag = "0"))
  ) |>
  filter(valid_name %in% rownames(trait)) |>
  dplyr::count(comm_code, valid_name, name = "count") |>
  tidyr::pivot_wider(names_from = valid_name, values_from = count, values_fill = 0) |>
  tibble::column_to_rownames(var = "comm_code") |>
  as.matrix()

# Validate
stopifnot(length(setdiff(colnames(comm), rownames(trait))) == 0)
message("Community matrix: ", nrow(comm), " sites × ", ncol(comm), " species")

saveRDS(comm, file.path(output_dir, "community_matrix.rds"))
saveRDS(trait, file.path(output_dir, "trait_matrix.rds"))
