# =============================================================================
# Script:  02e_assemble_diversity_table.R
# Purpose: Combine TD, SES-PD, SES-FD into the final diversity table and
#          compute µ-diversity.
#
# Input:   output/intermediate/sesPD.csv
#          output/intermediate/sesFD.csv
# Output:  data/processed/diversity_metrics_per_basin.csv
#
# Formula: µD = TD × ED, where ED = (zPD + zFD) / 2
#          zPD and zFD are cumulative-normal-rescaled SES values.
#          Reference: Tordoni et al. (2024) Glob. Ecol. Biogeogr.
# =============================================================================

library(tidyverse)
library(here)

output_dir <- here("output", "intermediate")
final_dir  <- here("data", "processed")

sesPD <- read.csv(file.path(output_dir, "sesPD.csv"))
sesFD <- read.csv(file.path(output_dir, "sesFD.csv"))

div_tbl <- left_join(sesFD, sesPD, by = c("SR","site")) |>
  mutate(
    # Cumulative normal rescaling (mean=0, SD=1)
    zPD  = as.numeric(pnorm(scale(sesPD), 0, 1)),
    zFD  = as.numeric(pnorm(scale(sesFD), 0, 1)),
    # Ecological distinctness
    ED   = rowMeans(cbind(zPD, zFD), na.rm = TRUE),
    # Taxonomic diversity (log species richness, standardised)
    logSR = log(SR),
    TD    = as.numeric(scale(logSR)),
    # µ-diversity
    muDiv = logSR * ED,
    # Basin code
    SBD_CD = as.numeric(str_remove(site, "Comm_"))
  ) |>
  dplyr::select(SBD_CD, SR, PD, sesPD, FD, sesFD, TD, zPD, zFD, ED, muDiv)

write.csv(div_tbl, file.path(final_dir, "diversity_metrics_per_basin.csv"), row.names = FALSE)
message("Final diversity table saved: ", nrow(div_tbl), " basins")
