# =============================================================================
# Script:  05b_figS4_distribution_maps.R
# Project: Amazon multidimensional biodiversity
# Purpose: Generate Supplementary Figure S4 — spatial distribution and
#          latitudinal gradient of each diversity metric (standardized only).
#
# Inputs:
#   data/processed/dataprepared_filtered.rds   (from 03b)
#   data/raw/SubDrainageShapefile/SubDrainageShapefile.shp
#   data/raw/env/habitat/riverine_network/ne_50m_rivers_lake_centerlines.shp
#   data/raw/env/anthropogenic/ne_cities_10.shp
#
# Outputs:
#   output/figures/si/FigS4_diversity_distribution_latitudinal.png
#
# Dependencies: tidyverse, sf, scales, paletteer, ggpubr, cowplot, here
# =============================================================================

library(tidyverse)
library(sf)
library(scales)
library(paletteer)
library(ggpubr)
library(cowplot)
library(here)

# --- Paths ---
input_dir   <- here("data", "processed")
basin_path  <- here("data", "raw", "SubDrainageShapefile", "SubDrainageShapefile.shp")
habitat_dir <- here("data", "raw", "env", "habitat")
anthro_dir  <- here("data", "raw", "env", "anthropogenic")
fig_si      <- here("output", "figures", "si")
dir.create(fig_si, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load basemap and data ---
message("Loading data for SI Fig. S4...")

basins <- st_read(basin_path, quiet = TRUE) |> st_make_valid()
basins <- basins[order(basins$SBD_CD), ]
rownames(basins) <- basins$SBD_CD

rivers50 <- st_read(
  file.path(habitat_dir, "riverine_network", "ne_50m_rivers_lake_centerlines",
            "ne_50m_rivers_lake_centerlines.shp"), quiet = TRUE
) |> st_make_valid()

dataprepared_filt <- readRDS(file.path(input_dir, "dataprepared_filtered.rds"))

amz_df <- basins |> left_join(dataprepared_filt, by = "SBD_CD")

# --- 2. Prepare latitudinal gradient data ---
# CENT_Y.x comes from the join (dataprepared has CENT_Y)
dbLat <- amz_df |>
  st_as_sf() |>
  as.data.frame() |>
  dplyr::select(SBD_CD, logSR, sesPD, sesFD, muDiv, CENT_Y.x) |>
  drop_na() |>
  group_by(SBD_CD) |>
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

river <- st_intersection(st_as_sf(rivers50), st_make_valid(amz_df))

# --- 3. Metric definitions ---
metrics <- c("logSR", "sesPD", "sesFD", "muDiv")

metric_palettes <- c(
  muDiv = "ggsci::red_material",
  logSR = "RColorBrewer::BuGn",
  sesPD = "MoMAColors::Avedon",
  sesFD = "MetBrewer::Benedictus"
)

metric_labels <- c(
  muDiv = "μ-diversity",
  logSR = "Taxonomic diversity",
  sesPD = "Phylogenetic diversity",
  sesFD = "Functional diversity"
)

# --- 4. Refactored map + gradient function ---
make_distribution_panel <- function(metric, amz_data, lat_data, river_data, city_data) {
  
  current_pal <- paletteer_d(metric_palettes[metric], n = 6)
  # Reverse gradient for sesFD and sesPD (darker = higher)
  pal_gradient <- if (metric %in% c("sesFD", "sesPD")) {
    c(current_pal[6:1])
  } else {
    c(current_pal[1:6])
  }
  response_label <- metric_labels[metric]
  
  # Spatial map
  dataAux <- amz_data |>
    st_as_sf() |>
    dplyr::select(SBD_CD, all_of(metric), CENT_Y.x, geometry)
  names(dataAux)[2] <- "varAux"
  
  map1 <- ggplot(data = dataAux, aes(fill = varAux)) +
    geom_sf(linewidth = 0.1) +
    geom_sf(data = river_data, inherit.aes = FALSE,
            colour = "dodgerblue", alpha = 0.5, size = 1) +
    do.call(scale_fill_gradientn, list(
      colors = pal_gradient,
      labels = comma_format(accuracy = 1),
      na.value = "grey70"
    )) +
    coord_sf(ylim = c(-20, 5), xlim = c(-79, -51), expand = TRUE) +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(vjust = 0.7, size = 16),
      legend.key.size = unit(0.4, "cm"),
      legend.direction = "horizontal",
      legend.box.just = "center",
      legend.margin = margin(t = -30, b = 15, unit = "pt"),
      panel.border = element_rect(colour = "black", fill = NA),
      plot.margin = unit(c(0.2, 0.2, 0.1, 0.2), "cm")
    ) +
    labs(fill = response_label) +
    guides(fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(5, "cm"), barheight = unit(0.2, "cm"),
      direction = "horizontal"
    ))
  
  # Latitudinal gradient
  latAux <- lat_data |>
    dplyr::select(SBD_CD, all_of(metric), CENT_Y.x) |>
    setNames(c("SBD_CD", "varAux", "LAT"))
  
  map_gradient <- ggplot(data = latAux, aes(x = LAT, y = varAux, color = varAux)) +
    geom_point(col = "gray70") +
    geom_smooth(aes(color = after_stat(y)), se = FALSE, method = "loess", linewidth = 3) +
    do.call(scale_color_gradientn, list(
      colors = pal_gradient,
      labels = comma_format(accuracy = 1),
      na.value = "grey70"
    )) +
    scale_x_continuous(breaks = seq(from = -55, to = 80, by = 20)) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.margin = unit(c(0.2, 0.2, 0.1, 0.2), "cm")
    ) +
    coord_flip() +
    labs(y = response_label, x = "Latitude (\u00b0)")
  
  # Combine: gradient (left) + map (right)
  ggarrange(map_gradient, map1, ncol = 2, widths = c(1, 2), align = "h")
}

# --- 5. Generate all four panels ---
message("Generating distribution panels...")

listaFigs <- lapply(metrics, function(m) {
  make_distribution_panel(m, amz_df, dbLat, river, city)
})

# --- 6. Combine into 2x2 figure and save ---
message("Assembling SI Fig. S4...")

figS4 <- ggarrange(
  plot_grid(listaFigs[[1]]),
  plot_grid(listaFigs[[2]]),
  plot_grid(listaFigs[[3]]),
  plot_grid(listaFigs[[4]]),
  ncol = 2, nrow = 2, labels = "AUTO"
)

ggsave(file.path(fig_si, "FigS4_diversity_distribution_latitudinal.png"),
       plot = figS4, width = 5000, height = 4000, units = "px", dpi = 400)

message("05b complete. SI Fig. S4 saved.")
