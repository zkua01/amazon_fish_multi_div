# =============================================================================
# Script:  01c_clean_occurrences.R
# Project: Amazon multidimensional biodiversity
# Purpose: Clean raw AmazonFish occurrence records, standardise species names,
#          and map to valid taxonomy using a reference lookup table.
#
# Input:   data/processed/amazon_fish_reftbl.csv          (reference lookup)
#          data/raw/Amazonfish-CompleteDatabase2022_v2.csv (public; see DOI below)
# Output:  data/processed/amazon_occurrences_clean.csv
#
# Note: The raw AmazonFish database is publicly available at:
#       https://doi.org/10.1038/s41597-020-0436-4
#       Download it before running this script.
#
# Dependencies: tidyverse, here
# =============================================================================

library(tidyverse)
library(here)

# --- Paths ---
input_dir  <- here("data", "processed")
raw_dir    <- here("data", "raw")
output_dir <- here("data", "processed")

# --- 1. Load reference lookup table ------------------------------------------
occ_fw_db <- read.csv(file.path(input_dir, "amazon_fish_reftbl.csv")) |>
  dplyr::select(occ_spp, valid_name) |>
  mutate(occ_spp = str_replace_all(occ_spp, " ", "_")) |>
  distinct()

# --- 2. Load raw occurrence database ------------------------------------------
occ_raw <- read.csv(
  file.path(raw_dir, "Amazonfish-CompleteDatabase2022_v2.csv"),
  sep = ";"
)

# --- 3. Helper: clean species name strings ------------------------------------
clean_names <- function(x) {
  x |>
    str_replace_all("\\.", " ") |>
    str_remove_all("[^A-Za-z_ ]") |>
    str_replace_all(" +", "_") |>
    str_trim() |>
    str_to_title()
}

# --- 4. Filter valid occurrences and clean names --------------------------------
occ_fw <- occ_raw |>
  filter(Occurrence.Status == "valid") |>
  mutate(
    Species   = clean_names(Referent.Species.Name),
    Longitude.X = as.numeric(str_replace(Longitude.X, ",", ".")),
    Latitude.Y  = as.numeric(str_replace(Latitude.Y, ",", "."))
  ) |>
  distinct()

rm(occ_raw)

# --- 5. Map to valid names via reference table ---------------------------------
occ_fw_clean <- occ_fw |>
  left_join(occ_fw_db, by = c("Species" = "occ_spp")) |>
  # Manual correction not in the lookup table
  mutate(
    valid_name = ifelse(
      valid_name == "Leporinus y-ophorus",
      "Leporinus yophorus",
      valid_name
    )
  )

# --- 6. Validate ----------------------------------------------------------------
n_total   <- nrow(occ_fw_clean)
n_unmapped <- sum(is.na(occ_fw_clean$valid_name))
message("Total valid occurrences: ", n_total)
message("Occurrences without valid_name mapping: ", n_unmapped)

# --- 7. Save ---------------------------------------------------------------------
write.csv(
  occ_fw_clean,
  file.path(output_dir, "amazon_occurrences_clean.csv"),
  row.names = FALSE
)
message("Cleaned occurrence data saved.")