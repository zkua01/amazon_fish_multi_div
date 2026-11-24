# =============================================================================
# Script:  06b_overlap_permutation.R
# Project: Amazon multidimensional biodiversity
# Purpose: Test whether the observed spatial decoupling between TD, PD, and FD
#          (quantified as low multi-facet overlap at the 17% threshold) is
#          significantly lower than expected by chance, using a rank-permutation
#          null model.
#
# Null model:
#   Hold the TD ranking fixed (as the reference axis). Independently permute
#   the PD and FD rankings across the 143 subdrainage basins. For each
#   permutation, compute the proportion of basins in the top 17% for all
#   four facets simultaneously (TD ∩ PD ∩ FD ∩ μD). Repeat 10,000 times.
#
#   Rationale for holding TD fixed:
#     TD (log species richness) is the metric most directly affected by
#     sampling effort. By anchoring the null to TD, we ask: "Given this
#     particular TD spatial pattern, how often would PD and FD coincidentally
#     align with it?" If the observed overlap is far below the null, the
#     decoupling cannot be explained by random co-occurrence of richness
#     and phylogenetic/functional structure.
#
#   Rationale for permuting PD and FD independently:
#     PD and FD are derived from the same species pool but capture orthogonal
#     information (evolutionary history vs. ecological traits). Independent
#     permutation tests whether their spatial concordance with TD exceeds
#     chance. A joint permutation (preserving PD–FD correlation) would test
#     a different (weaker) null; we use the conservative independent version.
#
#   μD is recomputed from the permuted PD and FD rankings (since μD = TD × ED,
#   where ED depends on zPD and zFD). This ensures the 4-way overlap is
#   internally consistent under the null.
#
# Inputs:
#   data/processed/diversity_metrics_per_basin.csv  (from 02e)
#
# Outputs:
#   output/intermediate/overlap_permutation_results.rds
#   Console summary (observed vs. null; p-value)
#
# Dependencies: tidyverse, here
# =============================================================================

library(tidyverse)
library(here)

# --- Paths ---
input_dir <- here("data", "processed")
inter_dir <- here("output", "intermediate")

# =============================================================================
# SECTION 1: LOAD DIVERSITY DATA
# =============================================================================

message("Loading diversity metrics...")
div_tbl <- read.csv(file.path(input_dir, "diversity_metrics_per_basin.csv")) |>
  mutate(logSR = log(SR))

# Verify we have the required columns
required_cols <- c("SBD_CD", "logSR", "sesPD", "sesFD", "muDiv")
stopifnot(all(required_cols %in% names(div_tbl)))

n_basins <- nrow(div_tbl)
message("Basins: ", n_basins)

# =============================================================================
# SECTION 2: COMPUTE OBSERVED OVERLAP
# =============================================================================

threshold <- 0.83  # Top 17% (Aichi Target 11 equivalent)

# Compute percent ranks
div_obs <- div_tbl |>
  mutate(
    rTD  = percent_rank(logSR),
    rPD  = percent_rank(sesPD),
    rFD  = percent_rank(sesFD),
    rMuD = percent_rank(muDiv)
  )

# Observed 4-way overlap
observed_overlap <- sum(
  div_obs$rTD >= threshold &
    div_obs$rPD >= threshold &
    div_obs$rFD >= threshold &
    div_obs$rMuD >= threshold,
  na.rm = TRUE
) / n_basins

message("Observed 4-way overlap at 17%: ", round(observed_overlap * 100, 2), "%")
message("  (", sum(div_obs$rTD >= threshold & div_obs$rPD >= threshold &
                     div_obs$rFD >= threshold & div_obs$rMuD >= threshold, na.rm = TRUE),
        " of ", n_basins, " basins)")

# Also compute pairwise overlaps for reference
pairwise_obs <- data.frame(
  pair = c("TD∩PD", "TD∩FD", "PD∩FD", "TD∩μD", "PD∩μD", "FD∩μD"),
  proportion = c(
    sum(div_obs$rTD >= threshold & div_obs$rPD >= threshold, na.rm = TRUE) / n_basins,
    sum(div_obs$rTD >= threshold & div_obs$rFD >= threshold, na.rm = TRUE) / n_basins,
    sum(div_obs$rPD >= threshold & div_obs$rFD >= threshold, na.rm = TRUE) / n_basins,
    sum(div_obs$rTD >= threshold & div_obs$rMuD >= threshold, na.rm = TRUE) / n_basins,
    sum(div_obs$rPD >= threshold & div_obs$rMuD >= threshold, na.rm = TRUE) / n_basins,
    sum(div_obs$rFD >= threshold & div_obs$rMuD >= threshold, na.rm = TRUE) / n_basins
  )
)
message("\nObserved pairwise overlaps:")
print(pairwise_obs)

# =============================================================================
# SECTION 3: PERMUTATION TEST
# =============================================================================

n_perm <- 10000
message("\nRunning ", n_perm, " permutations...")

set.seed(324)

# Pre-extract rankings as vectors for speed
rTD_vec  <- div_obs$rTD
rPD_vec  <- div_obs$rPD
rFD_vec  <- div_obs$rFD
logSR_vec <- div_obs$logSR

# Store results
null_overlaps_4way <- numeric(n_perm)
null_overlaps_TDPD <- numeric(n_perm)
null_overlaps_TDFD <- numeric(n_perm)
null_overlaps_PDFD <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  if (i %% 2000 == 0) message("  Permutation ", i, " / ", n_perm)
  
  # Permute PD and FD rankings independently
  perm_PD <- sample(rPD_vec)
  perm_FD <- sample(rFD_vec)
  
  # Recompute μD under the null:
  # ED = mean(zPD, zFD); μD = logSR × ED
  # Since we're permuting ranks, we use the permuted ranks directly
  # as proxies for zPD and zFD (monotonic transform preserves rank order)
  perm_ED <- (perm_PD + perm_FD) / 2
  perm_muD_rank <- percent_rank(logSR_vec * perm_ED)
  
  # Compute overlaps
  top_TD  <- rTD_vec >= threshold
  top_PD  <- perm_PD >= threshold
  top_FD  <- perm_FD >= threshold
  top_MuD <- perm_muD_rank >= threshold
  
  null_overlaps_4way[i] <- sum(top_TD & top_PD & top_FD & top_MuD) / n_basins
  null_overlaps_TDPD[i] <- sum(top_TD & top_PD) / n_basins
  null_overlaps_TDFD[i] <- sum(top_TD & top_FD) / n_basins
  null_overlaps_PDFD[i] <- sum(top_PD & top_FD) / n_basins
}

# =============================================================================
# SECTION 4: COMPUTE P-VALUES AND SUMMARY
# =============================================================================

p_4way <- mean(null_overlaps_4way >= observed_overlap)

message("\n=== Permutation test results ===")
message("Observed 4-way overlap: ", round(observed_overlap * 100, 2), "%")
message("Null 4-way overlap:     mean = ", round(mean(null_overlaps_4way) * 100, 2),
        "% ± ", round(sd(null_overlaps_4way) * 100, 2), "%")
message("Null range:             ", round(quantile(null_overlaps_4way, 0.025) * 100, 2),
        "% – ", round(quantile(null_overlaps_4way, 0.975) * 100, 2), "% (95% CI)")
message("p-value (one-sided):    ", p_4way)

if (p_4way < 0.001) {
  message("→ The observed decoupling is significant at p < 0.001.")
} else if (p_4way < 0.01) {
  message("→ The observed decoupling is significant at p < 0.01.")
} else {
  message("→ WARNING: The observed overlap is NOT significantly lower than random expectation.")
}

# Pairwise p-values
p_TDPD <- mean(null_overlaps_TDPD >= pairwise_obs$proportion[1])
p_TDFD <- mean(null_overlaps_TDFD >= pairwise_obs$proportion[2])
p_PDFD <- mean(null_overlaps_PDFD >= pairwise_obs$proportion[3])

message("\nPairwise permutation p-values:")
message("  TD∩PD: observed = ", round(pairwise_obs$proportion[1] * 100, 1),
        "%, null mean = ", round(mean(null_overlaps_TDPD) * 100, 1),
        "%, p = ", p_TDPD)
message("  TD∩FD: observed = ", round(pairwise_obs$proportion[2] * 100, 1),
        "%, null mean = ", round(mean(null_overlaps_TDFD) * 100, 1),
        "%, p = ", p_TDFD)
message("  PD∩FD: observed = ", round(pairwise_obs$proportion[3] * 100, 1),
        "%, null mean = ", round(mean(null_overlaps_PDFD) * 100, 1),
        "%, p = ", p_PDFD)

# =============================================================================
# SECTION 5: SAVE RESULTS
# =============================================================================

perm_results <- list(
  observed_overlap_4way = observed_overlap,
  null_overlaps_4way = null_overlaps_4way,
  p_value_4way = p_4way,
  pairwise_observed = pairwise_obs,
  null_overlaps_TDPD = null_overlaps_TDPD,
  null_overlaps_TDFD = null_overlaps_TDFD,
  null_overlaps_PDFD = null_overlaps_PDFD,
  p_TDPD = p_TDPD,
  p_TDFD = p_TDFD,
  p_PDFD = p_PDFD,
  n_perm = n_perm,
  threshold = threshold,
  n_basins = n_basins
)

saveRDS(perm_results, file.path(inter_dir, "overlap_permutation_results.rds"))
message("\nResults saved to output/intermediate/overlap_permutation_results.rds")
message("06b complete.")
