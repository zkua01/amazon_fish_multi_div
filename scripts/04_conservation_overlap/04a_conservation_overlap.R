# =============================================================================
# Script:  04a_conservation_overlap.R
# Project: Amazon multidimensional biodiversity
# Purpose: Quantify spatial overlap among diversity facets at multiple
#          conservation thresholds and assess alignment with existing
#          proactive conservation plans (Jézéquel et al. 2020).
#
# Outputs:
#   output/figures/main/Fig2_upset_overlap.png           (Fig. 2)
#   output/figures/extended/ExtendedDataFig_E3_overlap_maps.png  (Ext. Fig. E3)
#   output/figures/main/Fig6_conservation_overlap.png    (Fig. 6)
#
# Inputs:
#   data/processed/dataprepared_filtered.rds   (from 03b)
#   data/raw/SubDrainageShapefile/SubDrainageShapefile.shp
#   data/raw/env/habitat/riverine_network/ne_50m_rivers_lake_centerlines.shp
#   data/raw/env/anthropogenic/ne_cities_10.shp
#   data/processed/jezequel_overlap.csv        (conservation plan table)
#
# Dependencies: tidyverse, ComplexUpset, sf, eulerr, VennDiagram, ggplotify,
#               gridExtra, cowplot, patchwork, ggpubr, paletteer, here
# =============================================================================

library(tidyverse)
library(ComplexUpset)
library(sf)
library(scales)
library(patchwork)
library(cowplot)
library(VennDiagram)
library(eulerr)
library(ggplotify)
library(gridExtra)
library(ggpubr)
library(paletteer)
library(here)

# --- Paths ---
input_dir   <- here("data", "processed")
basin_path  <- here("data", "raw", "SubDrainageShapefile", "SubDrainageShapefile.shp")
habitat_dir <- here("data", "raw", "env", "habitat")
anthro_dir  <- here("data", "raw", "env", "anthropogenic")
fig_main    <- here("output", "figures", "main")
fig_ext     <- here("output", "figures", "extended")
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_ext,  showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SECTION 1: LOAD DATA
# =============================================================================

message("Loading basemap and ancillary data...")

# Basin shapefile
basins <- st_read(basin_path, quiet = TRUE)
basins <- basins[order(basins$SBD_CD), ]
rownames(basins) <- basins$SBD_CD

# Rivers (for basemap decoration)
rivers50 <- st_read(
  file.path(habitat_dir, "riverine_network", "ne_50m_rivers_lake_centerlines",
            "ne_50m_rivers_lake_centerlines.shp"),
  quiet = TRUE
) |> st_make_valid()

# Conservation plan table (Jézéquel et al. 2020)
conserve_tbl <- read.csv(file.path(input_dir, "jezequel_overlap.csv"))
proactive_2_cd <- conserve_tbl |>
  filter(plan == "pro_2") |> pull(SBD_CD) |>
  str_split(pattern = ",") |> unlist() |> str_squish()
proactive_3_cd <- conserve_tbl |>
  filter(plan == "pro_3") |> pull(SBD_CD) |>
  str_split(pattern = ",") |> unlist() |> str_squish()

# Diversity data (from 03b)
dataprepared_filt <- readRDS(file.path(input_dir, "dataprepared_filtered.rds"))

# Join basins with diversity data
amz_df <- basins |> left_join(dataprepared_filt, by = "SBD_CD")

# Intersect rivers and cities with basins (for map decoration)
river <- st_intersection(st_as_sf(rivers50), st_make_valid(amz_df))

message("Data loaded: ", nrow(amz_df), " basins")

# =============================================================================
# SECTION 2: COMPUTE PERCENT RANKS AND OVERLAP DATA
# =============================================================================

message("Computing percent ranks...")

Divdata <- amz_df |>
  dplyr::select(SBD_CD, logSR, sesPD, sesFD, muDiv) |>
  as.data.frame() |>
  dplyr::select(-geometry)

datapreparedTot <- Divdata |>
  mutate(
    rmuDiv = percent_rank(muDiv),
    rSR    = percent_rank(logSR),
    rPD    = percent_rank(sesPD),
    rFD    = percent_rank(sesFD),
    pro2   = ifelse(SBD_CD %in% proactive_2_cd, 1, 0),
    pro3   = ifelse(SBD_CD %in% proactive_3_cd, 1, 0)
  ) |>
  drop_na()

# Shared parameters
div_labels      <- c("μD", "FD", "PD", "TD")
thresholds      <- c(0.89, 0.83, 0.77, 0.50)
threshold_names <- c("Top 11%", "Top 17%", "Top 23%", "Top 50%")
degree_colors   <- c(`1` = "#ca0020", `2` = "#f4a582", `3` = "#92c5de", `4` = "#0571b0")

# Create boolean overlap data for each threshold
make_upset_data <- function(threshold) {
  datapreparedTot |>
    dplyr::select(-pro2, -pro3) |>
    drop_na() |>
    mutate(
      TD   = rSR    >= threshold,
      PD   = rPD    >= threshold,
      FD   = rFD    >= threshold,
      `μD` = rmuDiv >= threshold
    ) |>
    filter(TD | PD | FD | `μD`)
}

upsetdata_list <- lapply(thresholds, make_upset_data)
names(upsetdata_list) <- threshold_names

message("Overlap data computed for ", length(thresholds), " thresholds")

# =============================================================================
# SECTION 3: UPSET PLOTS (Figure 2)
# =============================================================================

message("Generating UpSet plots (Fig. 2)...")

# Note: query lists differ per threshold because different intersections
# are non-empty at each level. These are defined explicitly.

upset_t11_plot <- upset(upsetdata_list[["Top 11%"]], div_labels,
                        name = "Facet overlaps", wrap = TRUE, width_ratio = 0.2,
                        set_sizes = FALSE, sort_sets = FALSE,
                        base_annotations = list(
                          "Subbasin counts" = intersection_size(
                            text_colors = c(on_bar = "white", on_background = "black"),
                            text = aes(size = 2.5)
                          )
                        ),
                        queries = list(
                          upset_query(intersect = "μD", color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "TD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "PD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "FD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = c("μD", "TD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "PD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "FD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "FD", "PD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "FD", "TD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "FD", "PD", "TD"), color = degree_colors[4], fill = degree_colors[4])
                        ),
                        themes = upset_modify_themes(list(
                          "intersection_size" = theme(text = element_text(size = 0.5))
                        ))
) + ggtitle("Top 11%")

upset_t17_plot <- upset(upsetdata_list[["Top 17%"]], div_labels,
                        name = "Facet overlaps", wrap = TRUE, width_ratio = 0.2,
                        set_sizes = FALSE, sort_sets = FALSE,
                        base_annotations = list(
                          "Subbasin counts" = intersection_size(
                            text_colors = c(on_bar = "white", on_background = "black"),
                            text = aes(size = 2.5)
                          )
                        ),
                        queries = list(
                          upset_query(intersect = "μD", color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "TD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "PD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "FD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = c("μD", "PD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "TD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "FD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "FD", "PD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "TD", "PD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "FD", "TD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "FD", "PD", "TD"), color = degree_colors[4], fill = degree_colors[4])
                        ),
                        themes = upset_modify_themes(list(
                          "intersection_size" = theme(text = element_text(size = 0.5))
                        ))
) + ggtitle("Top 17%")

upset_t23_plot <- upset(upsetdata_list[["Top 23%"]], div_labels,
                        name = "Facet overlaps", wrap = TRUE, width_ratio = 0.2,
                        set_sizes = FALSE, sort_sets = FALSE,
                        base_annotations = list(
                          "Subbasin counts" = intersection_size(
                            text_colors = c(on_bar = "white", on_background = "black"),
                            text = aes(size = 2.5)
                          )
                        ),
                        queries = list(
                          upset_query(intersect = "μD", color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "TD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "PD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "FD",     color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = c("FD", "TD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("PD", "TD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "PD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "TD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "FD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "FD", "PD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "TD", "PD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "FD", "TD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "FD", "TD", "PD"), color = degree_colors[4], fill = degree_colors[4])
                        ),
                        themes = upset_modify_themes(list(
                          "intersection_size" = theme(text = element_text(size = 0.5))
                        ))
) + ggtitle("Top 23%")

upset_t50_plot <- upset(upsetdata_list[["Top 50%"]], div_labels,
                        name = "Facet overlaps", wrap = TRUE, width_ratio = 0.2,
                        set_sizes = FALSE, sort_sets = FALSE,
                        base_annotations = list(
                          "Subbasin counts" = intersection_size(
                            text_colors = c(on_bar = "white", on_background = "black"),
                            text = aes(size = 2.5)
                          )
                        ),
                        queries = list(
                          upset_query(intersect = "TD", color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "PD", color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = "FD", color = degree_colors[1], fill = degree_colors[1]),
                          upset_query(intersect = c("FD", "TD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("PD", "TD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("FD", "PD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "PD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "FD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "TD"), color = degree_colors[2], fill = degree_colors[2]),
                          upset_query(intersect = c("μD", "FD", "PD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "TD", "PD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "FD", "TD"), color = degree_colors[3], fill = degree_colors[3]),
                          upset_query(intersect = c("μD", "FD", "TD", "PD"), color = degree_colors[4], fill = degree_colors[4])
                        ),
                        themes = upset_modify_themes(list(
                          "intersection_size" = theme(text = element_text(size = 0.5))
                        ))
) + ggtitle("Top 50%")

# Combine and save (Figure 2)
upset_plot <- (upset_t11_plot + upset_t17_plot + upset_t23_plot + upset_t50_plot) +
  plot_layout(ncol = 2)

ggsave(file.path(fig_main, "Fig2_upset_overlap.png"),
       upset_plot, width = 5000, height = 4000, units = "px", dpi = 400)
message("Fig. 2 saved.")

# =============================================================================
# SECTION 4: DIVERSITY OVERLAP MAPS (Extended Figure E3)
# =============================================================================

message("Generating diversity overlap maps (Extended Fig. E3)...")

create_overlap_map <- function(data, title) {
  # Count how many facets each basin belongs to (1–4)
  data <- data |>
    mutate(n_facets = TD + PD + FD + `μD`) |>
    filter(n_facets > 0) |>
    mutate(n_facets = factor(n_facets, levels = 1:4))
  
  # Join with basemap geometry
  map_data <- basins |>
    left_join(data, by = "SBD_CD") |>
    filter(!is.na(n_facets))
  
  ggplot() +
    geom_sf(data = basins, fill = "white", color = "grey80", linewidth = 0.2) +
    geom_sf(data = map_data, aes(fill = n_facets), color = NA, alpha = 0.9) +
    geom_sf(data = river, colour = "dodgerblue", inherit.aes = FALSE,
            size = 1.5, alpha = 0.5) +
    scale_fill_manual(
      name = "Diversity overlap(s)",
      values = degree_colors,
      drop = FALSE
    ) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
}

map_t11 <- create_overlap_map(upsetdata_list[["Top 11%"]], "Top 11%")
map_t17 <- create_overlap_map(upsetdata_list[["Top 17%"]], "Top 17%")
map_t23 <- create_overlap_map(upsetdata_list[["Top 23%"]], "Top 23%")
map_t50 <- create_overlap_map(upsetdata_list[["Top 50%"]], "Top 50%")

div_overlap_map <- ggarrange(map_t11, map_t17, map_t23, map_t50,
                             common.legend = TRUE, legend = "bottom")

ggsave(file.path(fig_ext, "ExtendedDataFig_E3_overlap_maps.png"),
       div_overlap_map, width = 5000, height = 4000, units = "px", dpi = 400)
message("Extended Fig. E3 saved.")

# =============================================================================
# SECTION 5: CONSERVATION PLAN VENN DIAGRAMS (Figure 6 component)
# =============================================================================

message("Generating conservation plan Venn diagrams...")

# Colour palette for each diversity metric
col_pal <- c(
  md = paletteer_d("ggsci::red_material", n = 6)[6],
  td = paletteer_d("RColorBrewer::BuGn", n = 6)[6],
  pd = paletteer_d("MoMAColors::Avedon", n = 6)[1],
  fd = paletteer_d("MetBrewer::Benedictus", n = 6)[1]
)

# Refactored Venn function
make_venn <- function(metric_col, rank_col, threshold, metric_label, fill_colour) {
  venn_list <- list(
    metric = na.omit(datapreparedTot[datapreparedTot[[rank_col]] >= threshold, "SBD_CD"]),
    pro2   = na.omit(datapreparedTot[datapreparedTot$pro2 == 1, "SBD_CD"]),
    pro3   = na.omit(datapreparedTot[datapreparedTot$pro3 == 1, "SBD_CD"])
  )
  names(venn_list)[1] <- metric_label
  
  venndata <- VennDiagram::get.venn.partitions(venn_list) |>
    mutate(perc = round(..count.. / sum(..count..), 2)) |>
    dplyr::select(-..values..) |>
    setNames(c(metric_label, "Pro2", "Pro3", "set", "count", "perc"))
  
  venn_order <- c(
    paste0(metric_label, "&Pro2&Pro3"), "Pro2&Pro3",
    paste0(metric_label, "&Pro3"), "Pro3",
    paste0(metric_label, "&Pro2"), "Pro2", metric_label
  )
  
  venndiag <- venndata |>
    mutate(name = venn_order) |>
    pull(var = perc, name = name) |>
    euler(loss = "sum_absolute_region_error")
  
  as.ggplot(
    plot(venndiag,
         fills = list(fill = c(fill_colour, alpha("grey80", 0.5), alpha("grey15", 0.5))),
         quantities = list(cex = 0.5, type = "percent"),
         legend = NULL, bg = "transparent")
  ) + theme(plot.margin = unit(c(1, 1, 1, 1), "pt"))
}

vennplot_td <- make_venn("logSR", "rSR",    0.83, "TD",   col_pal[["td"]])
vennplot_pd <- make_venn("sesPD", "rPD",    0.83, "PD",   col_pal[["pd"]])
vennplot_fd <- make_venn("sesFD", "rFD",    0.83, "FD",   col_pal[["fd"]])
vennplot_mu <- make_venn("muDiv", "rmuDiv", 0.83, "μD", col_pal[["md"]])

message("Venn diagrams generated.")

# =============================================================================
# SECTION 6: CONSERVATION PLAN MAPS (Figure 6 component)
# =============================================================================

message("Generating conservation plan maps...")

# Prepare map data with top-17% indicators
compare_map_df <- datapreparedTot |>
  drop_na() |>
  mutate(
    muDiv17 = ifelse(rmuDiv >= 0.83, 1, 0),
    TD17    = ifelse(rSR    >= 0.83, 1, 0),
    PD17    = ifelse(rPD    >= 0.83, 1, 0),
    FD17    = ifelse(rFD    >= 0.83, 1, 0)
  ) |>
  left_join(basins, by = "SBD_CD") |>
  st_as_sf() |>
  st_make_valid()

# Refactored conservation map function
make_conserve_map <- function(indicator_col, metric_label, fill_colour, title) {
  ggplot(compare_map_df) +
    geom_sf(fill = "white", color = "grey80", linewidth = 0.2, alpha = 0.2) +
    geom_sf(data = filter(compare_map_df, .data[[indicator_col]] == 1),
            aes(fill = metric_label), color = NA) +
    geom_sf(data = filter(compare_map_df, pro2 == 1),
            aes(fill = "Proactive 2"), alpha = 0.5, color = NA) +
    geom_sf(data = filter(compare_map_df, pro3 == 1),
            aes(fill = "Proactive 3"), alpha = 0.5, color = NA) +
    coord_sf(expand = FALSE) +
    scale_fill_manual(
      name = "Conservation area",
      breaks = c(metric_label, "Proactive 2", "Proactive 3"),
      values = c(
        setNames(as.character(fill_colour), metric_label),
        "Proactive 2" = alpha("grey80", 0.5),
        "Proactive 3" = alpha("grey15", 0.5)
      )
    ) +
    labs(title = title) +
    theme_void() +
    theme(legend.position = "none", plot.margin = unit(c(0, 0, 0, 0), "pt"))
}

conserve_map_td <- make_conserve_map("TD17",    "TD",        col_pal[["td"]], "Taxonomic diversity")
conserve_map_pd <- make_conserve_map("PD17",    "PD",        col_pal[["pd"]], "Phylogenetic diversity")
conserve_map_fd <- make_conserve_map("FD17",    "FD",        col_pal[["fd"]], "Functional diversity")
conserve_map_mu <- make_conserve_map("muDiv17", "μD",   col_pal[["md"]], "μ-diversity")

message("Conservation maps generated.")

# =============================================================================
# SECTION 7: LEGEND AND COMBINED FIGURE 6
# =============================================================================

message("Assembling Fig. 6...")

# Shared legend
legend_plot <- ggplot(data.frame(variable = c("pro2", "pro3", "a", "b", "c", "d"))) +
  geom_point(aes(x = 1, y = 1, fill = factor(variable, levels = variable)),
             color = NA, show.legend = TRUE, size = 5, shape = 22) +
  scale_fill_manual(
    name = "Conservation priorities",
    values = c(
      "pro2" = alpha("grey80", 0.5),
      "pro3" = alpha("grey15", 0.5),
      "a" = as.character(col_pal["td"]),
      "b" = as.character(col_pal["pd"]),
      "c" = as.character(col_pal["fd"]),
      "d" = as.character(col_pal["md"])
    ),
    breaks = c("pro2", "pro3", "a", "b", "c", "d"),
    labels = c("Proactive 2", "Proactive 3",
               "Taxonomic diversity", "Phylogenetic diversity",
               "Functional diversity", "μ-diversity")
  ) +
  theme_void() +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.margin = margin(t = -10)
  )

# Overlay Venn diagrams onto maps
map_td <- ggdraw(conserve_map_td) + draw_plot(vennplot_td, x = 0.01, y = 0.01, width = 0.35, height = 0.35)
map_pd <- ggdraw(conserve_map_pd) + draw_plot(vennplot_pd, x = 0.01, y = 0.01, width = 0.35, height = 0.35)
map_fd <- ggdraw(conserve_map_fd) + draw_plot(vennplot_fd, x = 0.01, y = 0.01, width = 0.35, height = 0.35)
map_mu <- ggdraw(conserve_map_mu) + draw_plot(vennplot_mu, x = 0.01, y = 0.01, width = 0.35, height = 0.35)

# Arrange: 2x2 maps with Venn insets + legend at bottom
conserve_map <- arrangeGrob(map_td, map_pd, map_fd, map_mu,
                            ncol = 2, padding = unit(0, "line")) |>
  arrangeGrob(legend_plot, ncol = 1, heights = c(5, 1))

ggsave(file.path(fig_main, "Fig6_conservation_overlap.png"),
       grid.arrange(conserve_map), width = 5000, height = 4000, units = "px", dpi = 400)

message("04a complete.")
message("  - Fig. 2: ", file.path(fig_main, "Fig2_upset_overlap.png"))
message("  - Ext. Fig. E3: ", file.path(fig_ext, "ExtendedDataFig_E3_overlap_maps.png"))
message("  - Fig. 6: ", file.path(fig_main, "Fig6_conservation_overlap.png"))
