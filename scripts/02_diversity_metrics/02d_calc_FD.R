# =============================================================================
# Script:  02d_calc_FD.R
# Purpose: Calculate functional richness (TPD framework) and SES-FD.
#
# Input:   output/intermediate/community_matrix.rds
#          output/intermediate/trait_matrix.rds
#          output/intermediate/null_communities.rds
# Output:  output/intermediate/sesFD.csv
#
# Dependencies: TPD, Rcpp, parallel, custom C++ (src/fast_tpd_metrics_omp.cpp)
# =============================================================================

library(TPD)
library(tidyverse)
library(parallel)
library(Rcpp)
library(here)

output_dir <- here("output", "intermediate")
src_dir    <- here("src")

comm      <- readRDS(file.path(output_dir, "community_matrix.rds"))
trait     <- readRDS(file.path(output_dir, "trait_matrix.rds"))
null_comm <- readRDS(file.path(output_dir, "null_communities.rds"))

# Source custom functions
source(here("scripts", "utils", "R_fast_function.R"))
Rcpp::sourceCpp(file.path(src_dir, "fast_tpd_metrics_omp.cpp"))

# --- PCA on traits ---
PCAtraits <- princomp(scale(trait))
# Horn's parallel analysis retained 4 components
dim_used <- 4
traitUse <- as.matrix(PCAtraits$scores[, 1:dim_used])

# --- Community-level trait SD (for TPD kernels) ---
# Hscv estimation per community
Hpicom_raw <- matrix(NA, nrow(comm), dim_used)
for (i in seq_len(nrow(comm))) {
  spp <- colnames(comm)[which(comm[i, ] > 0)]
  if (length(spp) < 3) next
  traits_sub <- traitUse[spp, , drop = FALSE]
  hscv_est <- NULL
  for (k in dim_used:1) {
    hscv_est <- tryCatch(Hscv(traits_sub[, 1:k, drop = FALSE]), error = function(e) NULL)
    if (!is.null(hscv_est)) break
  }
  if (!is.null(hscv_est)) {
    n_used <- ncol(hscv_est)
    Hpicom_raw[i, 1:n_used] <- sqrt(diag(hscv_est))
  }
}

# Impute NAs with column means
col_means <- colMeans(Hpicom_raw, na.rm = TRUE)
Hpicom <- Hpicom_raw
for (j in seq_len(dim_used)) {
  Hpicom[is.na(Hpicom_raw[, j]), j] <- col_means[j]
}

# Uniform SD applied to all species (modelling simplification)
sdTraits <- matrix(rep(colMeans(na.omit(Hpicom)), nrow(traitUse)),
                   nrow = nrow(traitUse), byrow = TRUE, ncol = dim_used)
rownames(sdTraits) <- rownames(traitUse)

# --- Build sparse TPDs ---
fd_n_divisions <- 50
fd_tpd_threads <- max(1, detectCores() - 1)

fd_trait_ranges <- lapply(seq_len(dim_used), function(d)
  c(min(traitUse[, d] - 5 * sdTraits[, d]),
    max(traitUse[, d] + 5 * sdTraits[, d])))
fd_grid_list <- lapply(seq_len(dim_used), function(d)
  seq(fd_trait_ranges[[d]][1], fd_trait_ranges[[d]][2], length.out = fd_n_divisions))
fd_cell_volume <- prod(vapply(fd_grid_list, function(g) g[2] - g[1], numeric(1)))

sparse_TPDs_fd <- build_sparse_tpds_cpp(
  means = traitUse, sds = sdTraits, grid_list = fd_grid_list,
  alpha = 0.95, n_divisions = fd_n_divisions, n_threads = fd_tpd_threads
)

# --- Observed FD ---
sp_final <- intersect(colnames(comm), rownames(traitUse))
FDobs <- calc_fric_cpp(
  sparse_TPDs = sparse_TPDs_fd, abundance_mat = comm,
  cell_volume = fd_cell_volume, basin_names = rownames(comm)
)

# --- Null FD ---
calc_fd_null <- function(i) {
  simCom <- null_comm[[i]][, sp_final, drop = FALSE]
  calc_fric_cpp(
    sparse_TPDs = sparse_TPDs_fd, abundance_mat = simCom,
    cell_volume = fd_cell_volume, basin_names = rownames(simCom)
  )
}

FDrand <- sapply(seq_along(null_comm), calc_fd_null)

# --- SES-FD ---
fd_mean <- rowMeans(FDrand, na.rm = TRUE)
fd_sd   <- apply(FDrand, 1, sd, na.rm = TRUE)
fd_sd[fd_sd == 0 | is.na(fd_sd)] <- NA

sesFD <- data.frame(
  site  = rownames(comm),
  SR    = rowSums(comm > 0),
  FD    = FDobs,
  sesFD = (FDobs - fd_mean) / fd_sd
)

write.csv(sesFD, file.path(output_dir, "sesFD.csv"), row.names = FALSE)
message("SES-FD saved. Range: ", range(sesFD$sesFD, na.rm = TRUE))