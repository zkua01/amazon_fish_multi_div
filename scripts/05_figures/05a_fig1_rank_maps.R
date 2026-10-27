# =============================================================================
# Script:  05a_fig1_rank_maps.R
# Project: Amazon multidimensional biodiversity
# Purpose: Generate Figure 1 — subdrainage basins ranked by each diversity
#          metric (TD, PD, FD, μD) at six threshold categories.
#
# Inputs:
#   data/processed/dataprepared_filtered.rds   (from 03b)
#   data/raw/SubDrainageShapefile/SubDrainageShapefile.shp
#   data/raw/env/habitat/riverine_network/ne_50m_rivers_lake_centerlines.shp
#   data/raw/env/anthropogenic/ne_cities_10.shp
#
# Outputs:
#   output/figures/main/Fig1_diversity_rank_maps.png
#
# Dependencies: tidyverse, sf, paletteer, ggpubr, here
# =============================================================================

library(tidyverse)
library(sf)
library(paletteer)
library(ggpubr)
library(here)

# --- Paths ---
input_dir   <- here("data", "processed")
basin_path  <- here("data", "raw", "SubDrainageShapefile", "SubDrainageShapefile.shp")
habitat_dir <- here("data", "raw", "env", "habitat")
anthro_dir  <- here("data", "raw", "env", "anthropogenic")
fig_main    <- here("output", "figures", "main")
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load basemap and data ---
message("Loading data for Fig. 1...")

basins <- st_read(basin_path, quiet = TRUE) |> st_make_valid()
basins <- basins[order(basins$SBD_CD), ]
rownames(basins) <- basins$SBD_CD

rivers50 <- st_read(
  file.path(habitat_dir, "riverine_network", "ne_50m_rivers_lake_centerlines",
            "ne_50m_rivers_lake_centerlines.shp"), quiet = TRUE
) |> st_make_valid()

dataprepared_filt <- readRDS(file.path(input_dir, "dataprepared_filtered.rds"))

amz_df <- basins |> left_join(dataprepared_filt, by = "SBD_CD")

# --- 2. Compute percent ranks ---
amz_df_ranks <- amz_df |>
  mutate(
    rSR       = percent_rank(logSR),
    rPD       = percent_rank(sesPD),
    rFD       = percent_rank(sesFD),
    muDivRank = percent_rank(muDiv)
  )

river <- st_intersection(st_as_sf(rivers50), st_make_valid(amz_df_ranks))

# --- 3. Refactored rank map function ---
rank_levels <- c("Top 11%", "11%-17%", "17%-23%", "23%-50%", "50%-80%", "80%-100%")

make_rank_map <- function(data, rank_col, palette_name, title_text, direction = 1) {
  data_sf <- data |>
    st_as_sf() |>
    mutate(fct = case_when(
      .data[[rank_col]] >= 0.89 ~ "Top 11%",
      .data[[rank_col]] < 0.89 & .data[[rank_col]] >= 0.83 ~ "11%-17%",
      .data[[rank_col]] < 0.83 & .data[[rank_col]] >= 0.77 ~ "17%-23%",
      .data[[rank_col]] < 0.77 & .data[[rank_col]] >= 0.50 ~ "23%-50%",
      .data[[rank_col]] < 0.50 & .data[[rank_col]] >= 0.20 ~ "50%-80%",
      TRUE ~ "80%-100%"
    )) |>
    mutate(fct = factor(fct, levels = rank_levels))
  
  ggplot(data_sf, aes(fill = fct)) +
    geom_sf(data = data_sf |> filter(!is.na(.data[[rank_col]])),
            aes(fill = fct), linewidth = 0.05) +
    # Highlight basin with insufficient data (grey fill, red border)
    geom_sf(data = data_sf |> filter(is.na(.data[[rank_col]])),
            fill = "grey", color = "red", linewidth = 0.1) +
    geom_sf(data = river, colour = "dodgerblue", inherit.aes = FALSE,
            size = 1.5, alpha = 0.5) +
    coord_sf(ylim = c(-20, 5), xlim = c(-79, -51), expand = TRUE) +
    scale_fill_manual(values = paletteer_d(palette_name, n = 6, direction = direction),
                      na.translate = FALSE) +
    labs(y = "Latitude", x = "Longitude", fill = "Diversity ranks", title = title_text) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      legend.position = "none",
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.border = element_rect(fill = NA),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

# --- 4. Generate four rank maps ---
message("Generating rank maps...")

Fig1_SR <- make_rank_map(amz_df_ranks, "rSR",
                         "RColorBrewer::BuGn", "Taxonomic diversity", direction = -1)

Fig1_PD <- make_rank_map(amz_df_ranks, "rPD",
                         "MoMAColors::Avedon", "Phylogenetic diversity", direction = 1)

Fig1_FD <- make_rank_map(amz_df_ranks, "rFD",
                         "MetBrewer::Benedictus", "Functional diversity", direction = 1)

Fig1_MuDiv <- make_rank_map(amz_df_ranks, "muDivRank",
                            "ggsci::red_material", "\u03bc-diversity", direction = -1)

# --- 5. Combine and save ---
message("Assembling Fig. 1...")
fig1 <- ggarrange(Fig1_SR, Fig1_PD, Fig1_FD, Fig1_MuDiv,
                  legend = "right", labels = "AUTO")

ggsave(file.path(fig_main, "Fig1_diversity_rank_maps.png"),
       plot = fig1, width = 5000, height = 4000, units = "px", dpi = 400)

message("05a complete. Fig. 1 saved.")
