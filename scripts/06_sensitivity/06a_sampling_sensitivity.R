# =============================================================================
# Script:  06a_sampling_sensitivity.R
# Project: Amazon multidimensional biodiversity
# Purpose: Assess robustness of BRT results and conservation overlap patterns
#          to sampling effort inequality by excluding basins at the extremes
#          of sampling density.
#
# Design:
#   Four trimming levels are applied to the sampling effort predictor
#   (sample_area = unique sampling events per km²):
#
#
# Inputs:
#   data/processed/dataprepared_filtered.rds   (from 03b)
#   output/intermediate/brt_best_params.csv    (from 03c)
#   data/raw/SubDrainageShapefile/SubDrainageShapefile.shp
#
# Outputs:
#   output/intermediate/sensitivity_influence_comparison.csv
#   output/intermediate/sensitivity_overlap_comparison.csv
#   output/figures/extended/ExtendedDataFig_E8_sampling_sensitivity.png
#
# References:
#   Elith, J., Leathwick, J. R. & Hastie, T. (2008) A working guide to
#     boosted regression trees. J. Anim. Ecol. 77, 802–813.
#   Hastie, T., Tibshirani, R. & Friedman, J. (2009) The Elements of
#     Statistical Learning. Springer.
#   Pradervand, S. et al. (2014) Modelling the potential distribution of
#     European butterflies under climate change. Ecography 37, 1–11.
#   Tordoni, E. et al. (2024) µ-diversity: an integrative metric...
#     Glob. Ecol. Biogeogr. 33, e13845.
#   Jézéquel, D. et al. (2020) AmazonFish... Sci. Data 7, 163.
#
# Dependencies: gbm3, blockCV, sf, tidyverse, Rmisc, ComplexUpset, here
# =============================================================================

library(gbm3)
library(blockCV)
library(sf)
library(tidyverse)
library(Rmisc)
library(ComplexUpset)
library(here)

source(here("scripts", "utils", "brt_function.R"))

# --- Paths ---
input_dir <- here("data", "processed")
inter_dir <- here("output", "intermediate")
basin_path <- here("data", "raw", "SubDrainageShapefile", "SubDrainageShapefile.shp")
fig_si <- here("output", "figures", "si")
dir.create(fig_ext, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SECTION 1: LOAD DATA AND DEFINE TRIMMING
# =============================================================================

message("Loading data...")
dataprepared_filt <- readRDS(file.path(input_dir, "dataprepared_filtered.rds"))
best_params <- read.csv(file.path(inter_dir, "brt_best_params.csv"))
basins            <- st_read(basin_path, quiet = TRUE) |>
  st_transform(crs = 4326) |> arrange(SBD_CD)

responses  <- c("muDiv", "logSR", "sesFD", "sesPD")
resp_labels <- c(muDiv = "μ-diversity", logSR = "Taxonomic diversity",
                 sesFD = "Functional diversity", sesPD = "Phylogenetic diversity")

# Compute sampling effort quantiles
sq <- quantile(dataprepared_filt$sample_area, probs = c(0.05, 0.10, 0.15, 0.20, 0.80, 0.85, 0.90, 0.95), na.rm = TRUE)
message("Sampling effort quantiles (sites/km²):")
message("  5th: ", round(sq[1], 5))
message("  10th: ", round(sq[2], 5))
message("  15th: ", round(sq[3], 5))
message("  20th: ", round(sq[4], 5))
message("  80th: ", round(sq[5], 5))
message("  85th: ", round(sq[6], 5))
message("  90th: ", round(sq[7], 5))
message("  95th: ", round(sq[8], 5))

# Define trimming levels
trim_levels <- list(
  none = list(lo = -Inf, hi = Inf, label = "Full dataset (n=143)"),
  trim5 = list(lo = sq[1], hi = sq[8], label = "Remove 5% extremes (n=127)"),
  trim10 = list(lo = sq[2], hi = sq[7], label = "Remove 10% extremes (n=113)"),
  trim15 = list(lo = sq[3], hi = sq[6], label = "Remove 15% extremes (n=99)"),
  trim20 = list(lo = sq[4], hi = sq[5], label = "Remove 20% extremes (n=85)")
)

# =============================================================================
# SECTION 2: RE-RUN BRT FOR EACH TRIMMING LEVEL
# =============================================================================

# Shared fold creation function
create_spatial_folds <- function(data, cv, polygon) {
  data_sf <- st_as_sf(data, coords = c("CENT_X", "CENT_Y"), crs = 4326)
  cv_spatial(x = data_sf, k = cv, user_blocks = polygon,
             selection = "random", iteration = 100, progress = FALSE)$folds_list
}

# Predictor names (explicit, matching 03c)
predictor_names <- c(
  "AREA", "PC1_Water", "PC2_Temp", "PC4_Energy", "PC1_CurrPast",
  "PC1_Elev", "slr25m", "slr100m", "waterfalls", "NetworkDensity",
  "soilDiv", "LCDiv", "DistMouth_km", "sample_area"
)

run_brt_on_subset <- function(dp_subset, response, params_row, polygon) {
  "Fit BRT with fixed hyperparameters on a subset of basins.
  Returns relative influence (%) for each predictor."
  
  div <- dplyr::select(dp_subset, all_of(response))
  coord <- dplyr::select(dp_subset, CENT_X, CENT_Y)
  env <- dplyr::select(dp_subset, all_of(predictor_names)) |>
    mutate(across(AREA, log1p)) |>
    mutate(across(everything(), scale)) |>
    mutate(across(everything(), as.numeric))
  
  df <- bind_cols(div, env, coord)
  predictors <- names(env)
  
  # Fit model on full subset (no CV needed; we're comparing influence)
  model <- fit_brt_model(df, response, predictors, params_row)
  if (is.null(model)) return(NULL)
  
  # Relative influence
  rel_inf <- gbm3::relative_influence(model, num_trees = params_row$n_trees)
  total_inf <- sum(rel_inf[predictors], na.rm = TRUE)
  if (total_inf == 0) return(NULL)
  (rel_inf[predictors] / total_inf) * 100
}

# Execute for each trimming level and response
message("\n=== Running BRT sensitivity analysis ===")
influence_results <- list()

for (trim_name in names(trim_levels)) {
  trim <- trim_levels[[trim_name]]
  message("\n--- Trimming level: ", trim$label, " ---")
  
  # Subset data
  dp_sub <- dataprepared_filt |>
    filter(sample_area >= trim$lo, sample_area <= trim$hi) |>
    na.omit() |>
    left_join(dplyr::select(basins, SBD_CD, geometry), by = "SBD_CD")
  
  message("  Basins retained: ", nrow(dp_sub))
  
  polygon <- dplyr::select(dp_sub, SBD_CD, geometry) |>
    st_as_sf() |> st_make_valid()
  
  influence_results[[trim_name]] <- list()
  
  for (resp in responses) {
    params_row <- best_params |> filter(div == resp)
    message("  Fitting BRT for: ", resp)
    
    set.seed(324)
    inf <- run_brt_on_subset(dp_sub, resp, params_row, polygon)
    
    if (!is.null(inf)) {
      influence_results[[trim_name]][[resp]] <- data.frame(
        predictor = names(inf),
        influence = as.numeric(inf),
        response = resp,
        trim_level = trim_name,
        n_basins = nrow(dp_sub)
      )
    }
  }
}

# Combine all influence results
influence_all <- bind_rows(unlist(lapply(influence_results, function(x) {
  lapply(x, function(y) y)
}), recursive = FALSE)) |>
  mutate(trim_level = factor(trim_level, levels = c("none", "trim5", "trim10", "trim15", "trim20")))

write.csv(influence_all, file.path(inter_dir, "sensitivity_influence_comparison.csv"),
          row.names = FALSE)
message("\nInfluence comparison saved.")

# =============================================================================
# SECTION 3: COMPARE RELATIVE INFLUENCE (SPEARMAN RANK CORRELATION)
# =============================================================================

message("\n=== Comparing relative influence rankings ===")
rank_corrs <- data.frame()

for (resp in responses) {
  base <- influence_all |> filter(trim_level == "none", response == resp) |>
    arrange(desc(influence)) |> pull(predictor)
  
  for (trim_name in c("trim5", "trim10", "trim15", "trim20")) {
    comp <- influence_all |> filter(trim_level == trim_name, response == resp) |>
      arrange(desc(influence)) |> pull(predictor)
    
    rho <- cor(rank(base), rank(comp), method = "spearman")
    rank_corrs <- bind_rows(rank_corrs, data.frame(
      response = resp, trim_level = trim_name, spearman_rho = rho
    ))
  }
}

message("Spearman rank correlations (full vs. trimmed):")
print(rank_corrs)

# =============================================================================
# SECTION 4: COMPARE CONSERVATION OVERLAP (UpSet STATISTICS)
# =============================================================================

message("\n=== Comparing conservation overlap at 17% threshold ===")

compute_overlap_stats <- function(dp_sub, threshold = 0.83) {
  "Compute pairwise and 4-way overlap proportions at a given threshold."
  ranks <- dp_sub |>
    mutate(
      rTD  = percent_rank(logSR),
      rPD  = percent_rank(sesPD),
      rFD  = percent_rank(sesFD),
      rMuD = percent_rank(muDiv)
    ) |>
    filter(!is.na(rTD), !is.na(rPD), !is.na(rFD), !is.na(rMuD))
  
  top <- ranks |>
    mutate(
      TD   = rTD  >= threshold,
      PD   = rPD  >= threshold,
      FD   = rFD  >= threshold,
      MuD  = rMuD >= threshold
    )
  
  n <- nrow(top)
  overlap_4way <- sum(top$TD & top$PD & top$FD & top$MuD) / n
  overlap_TD_PD <- sum(top$TD & top$PD) / n
  overlap_TD_FD <- sum(top$TD & top$FD) / n
  overlap_PD_FD <- sum(top$PD & top$FD) / n
  
  data.frame(
    n_basins = n,
    overlap_4way = overlap_4way,
    overlap_TD_PD = overlap_TD_PD,
    overlap_TD_FD = overlap_TD_FD,
    overlap_PD_FD = overlap_PD_FD
  )
}

overlap_results <- list()
for (trim_name in names(trim_levels)) {
  trim <- trim_levels[[trim_name]]
  dp_sub <- dataprepared_filt |>
    filter(sample_area >= trim$lo, sample_area <= trim$hi) |>
    na.omit()
  overlap_results[[trim_name]] <- compute_overlap_stats(dp_sub) |>
    mutate(trim_level = trim_name, label = trim$label)
}

overlap_comparison <- bind_rows(overlap_results)
write.csv(overlap_comparison, file.path(inter_dir, "sensitivity_overlap_comparison.csv"),
          row.names = FALSE)

message("Overlap comparison (17% threshold):")
print(overlap_comparison |> select(label, n_basins, overlap_4way, overlap_TD_PD))

# =============================================================================
# SECTION 5: GENERATE EXTENDED DATA FIGURE E8
# =============================================================================

message("\n=== Generating Extended Data Fig. E8 ===")

# Panel A: Grouped bar plot of relative influence (top 6 predictors)
influence_plot <- influence_all |>
  group_by(response, trim_level) |>
  ungroup() |>
  left_join(
    tibble(abb = predictor_names,
           label = c("Area", "Water scarcity", "Solar energy input",
                     "Productivity seasonality", "Hist. climate warming",
                     "Terrain elevation", "SLR 25m", "SLR 100m",
                     "Waterfall counts", "Network density", "Soil diversity",
                     "Land cover diversity", "Distance to mouth (km)",
                     "Sampling effort")),
    by = c("predictor" = "abb")
  ) |>
  mutate(
    response_lab = factor(resp_labels[response], levels = c("Taxonomic diversity","Phylogenetic diversity","Functional diversity","μ-diversity")),
    trim_lab = factor(trim_level, levels = c("trim20","trim15","trim10", "trim5", "none"),
                      labels = c("Trim 20% (n=85)", "Trim 15% (n=99)", "Trim 10% (n=113)", "Trim 5% (n=127)", "Full (n=143)"))
  )

plot_sensitivity_influence <- function(data, response_name) {
  
  resp_data <- data |> filter(response == response_name)
  
  ggplot(resp_data, aes(x = reorder(label, influence), y = influence, colour = trim_lab,
                        fill = trim_lab)) +
    # Lollipop: line segment
    geom_linerange(aes(ymin = 0, ymax = influence),
                   position = position_dodge(width = 0.6),
                   linewidth = 1.2, alpha = 0.6) +
    # Dot at the tip
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    # 7.14% threshold line (100 / 14 predictors)
    geom_hline(yintercept = 7.14, linetype = "dashed", colour = "red",
               linewidth = 0.6) +
    coord_flip(ylim = c(0, 65)) +
    scale_colour_manual(name = "Dataset", values = c("#76DDE5","#05B5CC","#0084BD","#00459E","#062F65")) +
    scale_fill_manual(name = "Dataset", values = c("#76DDE5","#05B5CC","#0084BD","#00459E","#062F65")) +
    labs(
      title = resp_labels[response_name],
      x = NULL,
      y = "Relative influence (%)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9, colour = "black"),
      axis.text.x = element_text(size = 10, colour = "black"),
      axis.title = element_text(size = 11),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",   # hidden here; shown in combined figure
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
}

p_td <- plot_sensitivity_influence(influence_plot, "logSR")
p_pd <- plot_sensitivity_influence(influence_plot, "sesPD")
p_fd <- plot_sensitivity_influence(influence_plot, "sesFD")
p_md <- plot_sensitivity_influence(influence_plot, "muDiv")

fig_e8 <- ggpubr::ggarrange(p_td, p_pd, p_fd, p_md, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom", labels = "AUTO")

ggsave(file.path(fig_si, "SFig_E8_sampling_sensitivity.png"),
       plot = fig_e8, width = 5500, height = 5000, units = "px", dpi = 400)

message("Extended Data Fig. E8 saved.")

# influence_wide <- dplyr::select(influence_plot, -n_basins, -trim_lab, -response, -label) |>
#   pivot_wider(names_from = trim_level, values_from = influence) |>
#   mutate(delta_5  = trim5 - none,
#          delta_10 = trim10 - none) |>
#   pivot_longer(c(delta_5, delta_10), names_to = "trim", values_to = "delta") |>
#   mutate(trim_lab = factor(trim, levels = c("delta_5", "delta_10"), labels = c("Δ Full-5%", "Δ Full-10%")),
#          response_lab = factor(response_lab, levels = c("μ-diversity","Functional diversity","Phylogenetic diversity","Taxonomic diversity"), labels = c("μD","FD","PD","TD")))
# 
# heatmap_d5 <- ggplot(filter(influence_wide, trim_lab == "Δ Full-5%"), aes(x = predictor, y = response_lab, fill = delta)) +
#   geom_tile(colour = "white") +
#   geom_text(aes(label = round(delta, 1)), size = 4) +
#   scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
#                        midpoint = 0, name = "Δ influence (%)") +
#   labs(x = NULL, y = NULL, title = "Δ Full-5%") +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 11, colour = "black", face = "bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
#         legend.title = element_text(size = 11),
#         legend.text = element_text(size = 10),
#         plot.title = element_text(size = 14, face = "bold"))
# heatmap_d10 <- ggplot(filter(influence_wide, trim_lab == "Δ Full-10%"), aes(x = predictor, y = response_lab, fill = delta)) +
#   geom_tile(colour = "white") +
#   geom_text(aes(label = round(delta, 1)), size = 4) +
#   scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
#                        midpoint = 0, name = "Δ influence (%)") +
#   labs(x = NULL, y = NULL, title = "Δ Full-10%") +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 11, colour = "black", face = "bold"),
#         axis.text.x = element_text(angle = 45, hjust = 1, size = 10, colour = "black"),
#         legend.title = element_text(size = 11),
#         legend.text = element_text(size = 10),
#         plot.title = element_text(size = 14, face = "bold"))
# 
# heatmap_delta <- ggpubr::ggarrange(heatmap_d5, heatmap_d10, common.legend = TRUE, legend = "bottom")
# ggsave(file.path(fig_ext, "ExtendedDataFig_E8_sampling_sensitivity_heatmap.png"),
#        plot = heatmap_delta, width = 5500, height = 5000, units = "px", dpi = 400)

message("06a complete.")