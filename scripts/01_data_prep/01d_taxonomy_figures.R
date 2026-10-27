# =============================================================================
# Script:  01d_taxonomy_figures.R
# Project: Amazon multidimensional biodiversity
# Purpose: Generate Supplementary Figures S1 and S2 showing taxonomic
#          composition of the Amazonian freshwater fish assemblage.
#
# Input:   data/processed/trait_imputed_2492.csv
# Output:  output/figures/si/Figure_S1_taxonomy_summary.png
#          output/figures/si/Figure_S2_order_richness.png
#
# Dependencies: tidyverse, here
# =============================================================================

library(tidyverse)
library(here)

# --- Paths ---
input_dir  <- here("data", "processed")
figure_dir <- here("output", "figures", "si")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load imputed trait data (for taxonomy) ---------------------------------
trait_taxa <- read.csv(file.path(input_dir, "trait_imputed_2492.csv")) |>
  dplyr::select(Class, Order, Family, Genus, Species, Genus.species)

# --- 2. Figure S1: Taxonomic count by rank -------------------------------------
tax_summary <- trait_taxa |>
  summarise(
    Order   = n_distinct(Order),
    Family  = n_distinct(Family),
    Genus   = n_distinct(Genus),
    Species = n_distinct(Genus.species)
  ) |>
  pivot_longer(
    cols      = everything(),
    names_to  = "Taxonomy",
    values_to = "Count"
  ) |>
  mutate(Taxonomy = factor(Taxonomy, levels = c("Order", "Family", "Genus", "Species")))

fig_s1 <- ggplot(tax_summary, aes(x = Taxonomy, y = Count)) +
  geom_col(fill = "skyblue") +
  geom_text(
    aes(label = Count, vjust = ifelse(Count == max(Count), 1.1, -0.2)),
    size = 6
  ) +
  labs(x = "Taxonomic rank", y = "Count") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14)
  )

ggsave(
  file.path(figure_dir, "Figure_S1_taxonomy_summary.png"),
  plot = fig_s1, width = 5000, height = 4000, units = "px", dpi = 400
)

# --- 3. Figure S2: Richness per order ------------------------------------------
order_tbl <- trait_taxa |>
  group_by(Order) |>
  summarise(count = n(), .groups = "drop") |>
  mutate(Order = reorder(Order, count))

order_sum <- sum(order_tbl$count)

fig_s2 <- ggplot(order_tbl, aes(x = Order, y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(
    aes(label = count, hjust = ifelse(count == max(count), 1.1, -0.2)),
    size = 5
  ) +
  labs(x = "Order", y = "Observed richness") +
  scale_y_continuous(
    sec.axis = sec_axis(~ . / order_sum, name = "Relative richness")
  ) +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )

ggsave(
  file.path(figure_dir, "Figure_S2_order_richness.png"),
  plot = fig_s2, width = 5000, height = 4000, units = "px", dpi = 400
)

message("Supplementary Figures S1 and S2 saved.")