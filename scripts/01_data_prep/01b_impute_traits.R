# =============================================================================
# Script:  01b_impute_traits.R
# Project: Amazon multidimensional biodiversity
# Purpose: Impute missing functional trait values using phylogenetic
#          eigenvector regression (PVR) and random forest (missForest).
#
# Input:   data/processed/trait_miss_2492.csv
#          output/derived/phylogeny_amazon_2492.rds
# Output:  data/processed/trait_imputed_2492.csv
#
# Dependencies: PVR, missForest, tidyverse, here
# =============================================================================

library(PVR)
library(missForest)
library(FishPhyloMaker)
library(tidyverse)
library(here)

# --- Paths ---
input_dir  <- here("data", "processed")
tree_dir   <- here("output", "derived")
output_dir <- here("data", "processed")

# --- 1. Load trait data and phylogenetic tree --------------------------------
trait_raw <- read.csv(file.path(input_dir, "trait_miss_2492.csv"))
tree      <- readRDS(file.path(tree_dir, "phylogeny_amazon_2492.rds"))[["Phylogeny"]]

# --- 2. Add taxonomy from FishBase ------------------------------------------
tax_cache <- here::here("data", "processed", "fishbase_taxonomy_cache.csv")

if (file.exists(tax_cache)) {
  fb_taxa <- read.csv(tax_cache)
} else {
  fb_taxa_raw <- rfishbase::load_taxa(server = "fishbase")
  fb_taxa <- fb_taxa_raw |>
    select(Class, Order, Family, Genus) |>
    distinct()
  write.csv(fb_taxa, tax_cache, row.names = FALSE)
}

trait <- trait_raw |>
  mutate(Genus = word(Genus.species, 1),
         Species = word(Genus.species, 2)) |>
  left_join(fb_taxa, by = 'Genus') |>
  # fill in two remaining entries
  mutate(Family = case_when(Genus %in% c("Araichthys", "Loraxichthys","Lipopterichthys","Aposturisoma", "Etsaputu") ~ "Loricariidae",
                            Genus == "Melanosternarchus" ~ "Apteronotidae",
                            Genus == "Stegostenopos" ~ "Rhamphichthyidae",
                            Genus == "Piabucina" ~ "Lebiasinidae",
                            TRUE ~ Family),
         Order = case_when(Family == "Loricariidae" ~ "Siluriformes",
                           Family == "Apteronotidae" ~ "Gymnotiformes",
                           Family == "Rhamphichthyidae" ~ "Gymnotiformes",
                           Family == "Lebiasinidae" ~ "Characiformes",
                           TRUE ~ Order),
         Class = case_when(Order == "Siluriformes" ~ "Teleostei",
                           Order == "Gymnotiformes" ~ "Teleostei",
                           Order == "Characiformes" ~ "Teleostei",
                           TRUE ~ Class),
         MBl = ifelse(Genus.species == "Melanosternarchus amaru", 27.1, MBl),
         SpecCode = ifelse(Genus.species == "Melanosternarchus amaru", 71060, SpecCode)) |> 
  mutate(Genus.species = factor(Genus.species),
         SpecCode = factor(SpecCode),
         Genus = factor(Genus),
         Species = factor(Species),
         Class = factor(Class),
         Order = factor(Order),
         Family = factor(Family)) |>
  mutate(Genus.species = str_replace_all(Genus.species, " ", "_"),
         MBl = ifelse(MBl == -Inf, NA, MBl))

# Verify no missing taxonomy
stopifnot(!any(is.na(trait$Family) | is.na(trait$Order) | is.na(trait$Class)))

# --- 3. Filter to species present in the tree --------------------------------
trait_tree <- trait |>
  filter(Genus.species %in% tree$tip.label) |>
  mutate_if(is.numeric, log1p)  # log1p-transform for imputation

# --- 4. PVR decomposition ----------------------------------------------------
pvr_result <- PVRdecomp(tree, scale = TRUE)

# Retain top 10 eigenvectors (selected via scree plot inspection)
n_components <- 10
eigenvec <- data.frame(
  Genus.species = tree$tip.label,
  pvr_result@Eigen$vectors[, 1:n_components]
) |> distinct()

trait_miss <- full_join(trait_tree, eigenvec, by = "Genus.species") |>
  dplyr::select(-SpecCode, -Genus, -Species, -Family, -Class, -Order) |>
  mutate(Genus.species = factor(Genus.species))

# --- 5. Impute with missForest -----------------------------------------------
set.seed(324)
imputation <- missForest(trait_miss)

ximp <- imputation$ximp |>
  dplyr::select(-c(paste0("c", 1:n_components))) |>
  mutate_if(is.numeric, expm1)  # reverse log1p

# --- 6. Merge back and validate -----------------------------------------------
trait_fill <- trait |>
  dplyr::select(Class, Order, Family, Genus, Species, Genus.species, SpecCode, MBl:CPt) |>
  rows_patch(ximp, by = "Genus.species", unmatched = "ignore")

n_before <- nrow(trait_fill)
trait_fill <- trait_fill |> filter(!is.na(REs))
message("Species retained: ", nrow(trait_fill),
        " / ", n_before)

# --- 7. Save -------------------------------------------------------------------
write.csv(trait_fill,
          file.path(output_dir, "trait_imputed_2492.csv"),
          row.names = FALSE)
message("Imputed trait data saved.")