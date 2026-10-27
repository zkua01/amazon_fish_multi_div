# =============================================================================
# Script:  03b_variable_selection.R
# Project: Amazon multidimensional biodiversity
# Purpose: Select the final set of environmental predictors for BRT models
#          via Pearson correlation filtering (|r| >= 0.7) and VIF analysis.
#          Also computes spatial-corrected correlations between response
#          variables and generates Extended Data Figures E1 and E2.
#
# Inputs:
#   data/processed/dataprepared.rds                (from 03a)
#   data/processed/environmental_predictors.csv    (from 03a)
#
# Outputs:
#   data/processed/dataprepared_filtered.rds
#   output/figures/extended/ExtendedDataFig_E1_predictor_correlation_raw.png
#   output/figures/extended/ExtendedDataFig_E1_predictor_correlation_filtered.png
#   output/figures/extended/ExtendedDataFig_E2_response_correlation.png
#
# Variable selection summary (see Supplementary Table S2):
#   Dropped by |r| >= 0.7:  PC1_Temp, PC1_Energy, PC2_CurrPast, TributaryPosition
#   Dropped by VIF:         PC2_Energy, PC3_Energy, PC2_Water, WaterColor, PebasLake
#   Total predictors retained: 14
#
# Dependencies: ggcorrplot2, psych, tidyverse, SpatialPack, car, GGally, here
# =============================================================================

library(ggcorrplot2)
library(psych)
library(tidyverse)
library(SpatialPack)
library(GGally)
library(car)
library(here)

# --- Paths ---
input_dir  <- here("data", "processed")
figure_dir <- here("output", "figures", "extended")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load data ---
message("Loading prepared data...")
dataprepared <- readRDS(file.path(input_dir, "dataprepared.rds"))
envdata_full <- read.csv(file.path(input_dir, "environmental_predictors.csv"))

# --- 2. Extract response variables ---
# SR:muDiv selects all diversity metric columns (SR, FD, sesFD, PD, sesPD,
# zPD, zFD, ED, logSR, TD, muDiv)
response_df <- dplyr::select(dataprepared, c(SR:muDiv)) |>
  mutate(logSR = log(SR)) |>
  relocate(logSR, .after = ED) |>
  mutate(across(c(SR:muDiv), as.numeric)) |>
  drop_na()

# --- 3. Extract predictor variables ---
# Remove SBD_CD, CENT_X, CENT_Y (first 3 columns of envdata_full)
# Keep sample_area as a predictor
predictor_df_raw <- envdata_full[-c(1:3)]

# Log-transform area, then standardise continuous predictors
predictor_df <- predictor_df_raw |>
  mutate(LogArea = log1p(AREA)) |>
  dplyr::select(-AREA) |>
  mutate(across(!PebasLake & !WaterColor, scale),
         across(!PebasLake & !WaterColor, as.numeric))

# Classify variable types
continuous_vars  <- colnames(dplyr::select(predictor_df, -WaterColor, -PebasLake))
binary_vars      <- "PebasLake"
categorical_vars <- "WaterColor"

message("Predictors loaded: ", ncol(predictor_df), " total (",
        length(continuous_vars), " continuous, ",
        length(binary_vars), " binary, ",
        length(categorical_vars), " categorical)")

# =============================================================================
# SECTION A: PREDICTOR VARIABLE SELECTION
# =============================================================================

# --- 4. Pearson correlation filtering (|r| >= 0.7) ---
message("Checking predictor correlations...")
cor_tbl_pred <- corr.test(
  dplyr::select(predictor_df, -WaterColor, -PebasLake),
  method = "pearson"
)
cor_pred_r <- cor_tbl_pred$r
cor_pred_p <- cor_tbl_pred$p

# Identify highly correlated pairs
high_cor_pred_r <- as.data.frame(cor_pred_r) |>
  rownames_to_column(var = "var1") |>
  pivot_longer(cols = -var1, names_to = "var2", values_to = "r") |>
  filter(abs(r) >= 0.7, var1 != var2) |>
  rowwise() |>
  mutate(pair = paste(sort(c(var1, var2)), collapse = " x ")) |>
  ungroup() |>
  distinct(pair, .keep_all = TRUE) |>
  dplyr::select(-pair) |>
  arrange(desc(abs(r)))

message("Highly correlated pairs (|r| >= 0.7): ", nrow(high_cor_pred_r))
print(high_cor_pred_r)

# Extended Data Figure E1 (raw predictor correlation)
pred_corrplot <- ggcorrplot.mixed(
  cor_pred_r, upper = "circle", lower = "number",
  p.mat = cor_pred_p, insig = "label_sig",
  sig.lvl = c(0.05, 0.01, 0.001)
)
ggsave(file.path(figure_dir, "ExtendedDataFig_E1_predictor_correlation_raw.png"),
       plot = pred_corrplot, width = 5000, height = 4000, units = "px", dpi = 400)

# --- 5. Drop variables by correlation ---
# Decision based on domain knowledge and redundancy (see Supplementary Table S2)
drop_by_correlation <- c("PC1_Temp", "PC1_Energy", "PC2_CurrPast", "TributaryPosition")
message("Dropped by correlation: ", paste(drop_by_correlation, collapse = ", "))

filt_pred_df <- predictor_df |>
  dplyr::select(-all_of(drop_by_correlation))

# --- 6. VIF filtering ---
message("Checking VIF after correlation-based removal...")

# Build modelling dataframe for VIF
# Note: uses TD (scaled logSR) and muDiv, matching the original workflow
mod_data_filt <- dplyr::select(response_df, muDiv, sesFD, sesPD, TD) |>
  rownames_to_column(var = "SBD_CD") |>
  left_join(rownames_to_column(filt_pred_df, var = "SBD_CD"), by = "SBD_CD") |>
  na.omit()

filt_mod_md <- lm(muDiv ~ ., data = dplyr::select(mod_data_filt, -SBD_CD, -sesFD, -sesPD, -TD))
filt_mod_md_vif <- car::vif(filt_mod_md) |>
  as.data.frame() |>
  arrange(desc(GVIF))

message("VIF scores (top):")
print(head(filt_mod_md_vif, 10))

# Diagnostic: check which variables drive DistMouth_km and PC2_Temp collinearity
distmouth_mod <- lm(DistMouth_km ~ ., data = filt_pred_df)
message("DistMouth_km significant predictors:")
print(summary(distmouth_mod)$coefficients |>
        as.data.frame() |>
        arrange(-abs(`t value`)) |>
        filter(`Pr(>|t|)` < 0.05))

pc2temp_mod <- lm(PC2_Temp ~ ., data = filt_pred_df)
message("PC2_Temp significant predictors:")
print(summary(pc2temp_mod)$coefficients |>
        as.data.frame() |>
        arrange(-abs(`t value`)) |>
        filter(`Pr(>|t|)` < 0.05))

# Decision: drop variables with high VIF and/or categorical collinearity
drop_by_vif <- c("PC2_Energy", "PC3_Energy", "PC2_Water", "WaterColor", "PebasLake")
message("Dropped by VIF: ", paste(drop_by_vif, collapse = ", "))

# Verify VIF after second round of drops
filt_mod_md_2 <- lm(muDiv ~ ., data = dplyr::select(mod_data_filt,
                                                    -SBD_CD, -sesFD, -sesPD, -TD,
                                                    -PC2_Energy, -PC3_Energy, -PC2_Water, -WaterColor, -PebasLake))
filt_mod_md_2_vif <- car::vif(filt_mod_md_2) |> sort(decreasing = TRUE)
message("VIF after second round (should all be acceptable):")
print(filt_mod_md_2_vif)

# --- 7. Apply final variable selection ---
filt_pred_df_final <- filt_pred_df |>
  dplyr::select(-all_of(drop_by_vif))

message("Retained predictors (", ncol(filt_pred_df_final), "): ",
        paste(colnames(filt_pred_df_final), collapse = ", "))

# Extended Data Figure E1 (filtered predictor correlation)
cor_tbl_pred_filt <- corr.test(filt_pred_df_final, method = "pearson")
pred_filt_corrplot <- ggcorrplot.mixed(
  cor_tbl_pred_filt$r, upper = "circle", lower = "number",
  p.mat = cor_tbl_pred_filt$p, insig = "label_sig",
  sig.lvl = c(0.05, 0.01, 0.001)
)
ggsave(file.path(figure_dir, "ExtendedDataFig_E1_predictor_correlation_filtered.png"),
       plot = pred_filt_corrplot, width = 5000, height = 4000, units = "px", dpi = 400)

# =============================================================================
# SECTION B: RESPONSE CORRELATION (SPATIAL-CORRECTED)
# =============================================================================

# --- 8. Spatial-corrected correlation across all 7 response variables ---
message("Computing spatial-corrected response correlations...")
set.seed(324)
coords <- cbind(dataprepared$CENT_X, dataprepared$CENT_Y)

response_facet <- list(
  SR    = response_df$SR,
  PD    = response_df$PD,
  FD    = response_df$FD,
  logSR = response_df$logSR,
  sesPD = response_df$sesPD,
  sesFD = response_df$sesFD,
  muDiv = response_df$muDiv
)

res_pairs <- combn(names(response_facet), 2, simplify = FALSE)
corr_ttest_results <- lapply(res_pairs, function(p) {
  x <- response_facet[[p[[1]]]]
  y <- response_facet[[p[[2]]]]
  res <- modified.ttest(x, y, coords)
  data.frame(var1 = p[[1]], var2 = p[[2]],
             r = res$corr, p.value = res$p.value, ess = res$ESS)
}) |> bind_rows()

# Pearson correlations (for reference)
raw_cor_tbl_res <- corr.test(dplyr::select(response_df, SR, PD, FD), method = "pearson")
message("Raw metric correlations (SR, PD, FD):")
print(raw_cor_tbl_res, short = FALSE)

cor_tbl_res <- corr.test(dplyr::select(response_df, logSR, sesFD, sesPD, muDiv), method = "pearson")
message("Standardized metric correlations (logSR, sesFD, sesPD, muDiv):")
print(cor_tbl_res, short = FALSE)

# --- 9. Build correlation matrices for the 4 standardized metrics ---
# Define the 4 standardized metrics and their display labels (in desired order)
std_metrics   <- c("logSR", "sesFD", "sesPD", "muDiv")
display_order <- c("TD", "FD", "PD", "μD")  # μD

# Build symmetric pairwise results
corr_symmetric <- bind_rows(
  corr_ttest_results,
  rename(corr_ttest_results, var1 = var2, var2 = var1)
)

# --- Correlation matrix (r values) ---
cor_res_r <- corr_symmetric |>
  dplyr::select(-p.value, -ess) |>
  pivot_wider(names_from = var2, values_from = r, values_fill = 1) |>
  # (a) Filter rows to only the 4 standardized metrics
  filter(var1 %in% std_metrics) |>
  # (b) Select and rename columns
  dplyr::select(var1, TD = logSR, FD = sesFD, PD = sesPD, `μD` = muDiv) |>
  # Rename rows to display labels
  mutate(var1 = case_when(
    var1 == "logSR" ~ "TD",
    var1 == "sesFD" ~ "FD",
    var1 == "sesPD" ~ "PD",
    var1 == "muDiv" ~ "μD"
  )) |>
  # (c) Force row order to match column order
  mutate(var1 = factor(var1, levels = display_order)) |>
  arrange(var1) |>
  column_to_rownames("var1") |>
  as.matrix()

# --- P-value matrix ---
cor_res_p <- corr_symmetric |>
  dplyr::select(-r, -ess) |>
  pivot_wider(names_from = var2, values_from = p.value, values_fill = 0) |>
  filter(var1 %in% std_metrics) |>
  dplyr::select(var1, TD = logSR, FD = sesFD, PD = sesPD, `μD` = muDiv) |>
  mutate(var1 = case_when(
    var1 == "logSR" ~ "TD",
    var1 == "sesFD" ~ "FD",
    var1 == "sesPD" ~ "PD",
    var1 == "muDiv" ~ "μD"
  )) |>
  mutate(var1 = factor(var1, levels = display_order)) |>
  arrange(var1) |>
  column_to_rownames("var1") |>
  as.matrix()

# --- 10. Extended Data Figure E2: Response correlation ---
message("Generating Extended Data Fig. E2...")
res_corrplot <- ggcorrplot.mixed(
  cor_res_r, upper = "circle", lower = "number",
  p.mat = cor_res_p, insig = "label_sig",
  sig.lvl = c(0.05, 0.01, 0.001), pch.cex = 6
)
ggsave(file.path(figure_dir, "ExtendedDataFig_E2_response_correlation.png"),
       plot = res_corrplot, width = 4000, height = 3000, units = "px", dpi = 400)

# =============================================================================
# SECTION C: SAVE FILTERED DATASET
# =============================================================================

# --- 11. Assemble and save ---
all_dropped <- c(drop_by_correlation, drop_by_vif)
dataprepared_filt <- dataprepared |>
  mutate(logSR = log(SR)) |>
  relocate(logSR, .after = ED) |>
  dplyr::select(-all_of(all_dropped))

saveRDS(dataprepared_filt, file.path(input_dir, "dataprepared_filtered.rds"))

message("03b complete.")
message("  Filtered dataset: ", nrow(dataprepared_filt), " basins x ",
        ncol(dataprepared_filt), " columns")
message("  Dropped predictors (", length(all_dropped), "): ",
        paste(all_dropped, collapse = ", "))
message("  Output: data/processed/dataprepared_filtered.rds")