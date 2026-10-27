# =============================================================================
# Script:  02b_generate_null_communities.R
# Purpose: Generate 999 Curveball null communities (fixed row & column sums).
#
# Input:   output/intermediate/community_matrix.rds
# Output:  output/intermediate/null_communities.rds
#
# Reference: Strona et al. (2014) Nat. Commun. 5:4114
# =============================================================================

library(tidyverse)
library(parallel)
library(here)

source(here("scripts", "utils", "curveball.R"))

output_dir <- here("output", "intermediate")
comm <- readRDS(file.path(output_dir, "community_matrix.rds"))

set.seed(324)
perm   <- 999
ncore  <- max(1, detectCores() - 1)
cl     <- makeCluster(ncore)
clusterSetRNGStream(cl, 324)
clusterExport(cl, varlist = c("comm", "curve_ball"), envir = environment())

message("Generating ", perm, " null communities on ", ncore, " cores...")
null_comm <- parLapply(cl, seq_len(perm), function(i) curve_ball(comm))
stopCluster(cl)

saveRDS(null_comm, file.path(output_dir, "null_communities.rds"))
message("Done. Null communities saved.")
