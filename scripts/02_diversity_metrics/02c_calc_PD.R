# =============================================================================
# Script:  02c_calc_PD.R
# Purpose: Calculate Faith's PD and SES-PD using Curveball null model.
#
# Input:   output/intermediate/community_matrix.rds
#          output/intermediate/null_communities.rds
#          output/derived/phylogeny_amazon_2492.rds
# Output:  output/intermediate/sesPD.csv
# =============================================================================

library(tidyverse)
library(picante)
library(here)

output_dir <- here("output", "intermediate")
tree_dir   <- here("output", "derived")

comm      <- readRDS(file.path(output_dir, "community_matrix.rds"))
null_comm <- readRDS(file.path(output_dir, "null_communities.rds"))
tree      <- readRDS(file.path(tree_dir, "phylogeny_amazon_2492.rds"))[["Phylogeny"]]

# Clean tree tips
source(here("scripts", "utils", "clean_name.R"))
tree$tip.label <- clean_name(tree$tip.label)

# Observed PD (include.root = TRUE for consistency with null model)
PDobs <- pd(samp = comm, tree = tree, include.root = TRUE)$PD

# Null PD
PDrand <- matrix(NA, nrow = nrow(comm), ncol = length(null_comm))
rownames(PDrand) <- rownames(comm)
for (sim in seq_along(null_comm)) {
  simCom <- as.matrix(null_comm[[sim]])
  PDrand[, sim] <- pd(samp = simCom, tree = tree, include.root = TRUE)$PD
}

# SES
pd_null_mean <- rowMeans(PDrand, na.rm = TRUE)
pd_null_sd   <- apply(PDrand, 1, sd, na.rm = TRUE)
sesPD_vals   <- (PDobs - pd_null_mean) / pd_null_sd

sesPD <- data.frame(
  site  = rownames(comm),
  SR    = specnumber(comm),
  PD    = PDobs,
  sesPD = round(sesPD_vals, 2)
)

write.csv(sesPD, file.path(output_dir, "sesPD.csv"), row.names = FALSE)
message("SES-PD saved. Range: ", range(sesPD$sesPD))