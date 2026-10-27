# Code repository: Decoupling biodiversity dimensions reveals key biodiversity areas in the Amazon

**Authors:** Zi Xun KUA, Zhixuan SONG, Chengfan LI, Sébastien BROSSE, Pablo TEDESCO, Meng XU, Jun XU, Guohuan SU

---

## Overview

This repository contains the complete analysis pipeline for the manuscript *"Decoupling biodiversity dimensions reveals key biodiversity areas in the Amazon"*. The study demonstrates that taxonomic (TD), phylogenetic (PD), and functional (FD) diversity of Amazonian freshwater fishes are strongly spatially decoupled and governed by distinct assembly processes, with direct consequences for conservation prioritization.

---

## Repository Structure

```
amazon_fish_multi_div/
│
├── data/ # Gitignored; archived on Zenodo
│   ├── raw/
│   │   ├── SubDrainageShapefile/ # AmazonFish subdrainage basin polygons
│   │   └── env/
│   │       ├── climate/ # WorldClim, MODIS NPP, PET/AET, runoff
│   │       ├── habitat/ # Elevation, river network, runoff
│   │       └── historical/ # Paleoclimate rasters
│   ├── processed/ # Derived datasets (diversity metrics, env predictors, conservation plan)
│   └── spatial/
│
├── output/ # Gitignored; archive on Zenodo
│   ├── derived/ # Phylogenetic tree
│   ├── intermediate/ # PCA scores, BRT results, null communities
│   └── figures/
|       ├── main/ # Figs. 1-6
│       ├── extended/ # Extended Data Figures E1–E7
│       └── si/ # Supplementary Figures S1–S4
│
├── scripts/
│   ├── 01_data_prep/ # Biological and environmental data preparation
│   ├── 02_diversity_metrics/ # PD, FD, SES calculations, µ-diversity
│   ├── 03_brt_models/ # BRT fitting, bootstrap, validation, figures
│   ├── 04_conservation_overlap/ # UpSet plots, conservation plan alignment
│   ├── 05_figures/ # Fig. 1 and SI Fig. S4
│   └── utils/ # Shared utility functions
│
├── src/ # C++ source code for TPD calculations
└── README.md
```

---

## Script Descriptions and execution order

Scripts are designed be run sequentially within each directory. The recommended sequence is also ordered through the directory labels (01 → 05).

### `01_data_prep/` — Biological and environmental data preparation

| Script | Purpose | Key output |
|--------|---------|------------|
| `01a_make_phylo_tree.R` | Construct phylogenetic tree via FishPhyloMaker | `output/derived/phylogeny_amazon_2492.rds` |
| `01b_impute_traits.R` | Impute missing traits via PVR + missForest | `data/processed/trait_imputed_2492.csv` |
| `01c_clean_occurrences.R` | Clean AmazonFish occurrence records, map valid names | `data/processed/amazon_occurrences_clean.csv` |
| `01d_taxonomy_figures.R` | Generate SI Figs. S1–S2 (taxonomic composition) | `output/figures/si/Figure_S1_*.png`, `Figure_S2_*.png` |
| `01e_extract_climate_rasters.R` | Extract WorldClim, NPP, runoff; PCA (Water, Temp, Energy) | `output/intermediate/pca_water.csv`, `pca_temp.csv`, `pca_energy.csv` |
| `01f_process_historical_climate.R` | Extract paleoclimate rasters; climate anomaly PCA | `output/intermediate/pca_currpast.csv` |
| `01g_process_habitat_position.R` | Elevation, SLR, slope, network density, soil/LC diversity, position, sampling effort | `output/intermediate/pca_elev.csv`, `habitat_position_vars.csv` |

### `02_diversity_metrics/` — Diversity metric calculation

| Script | Purpose | Key output |
|--------|---------|------------|
| `02a_build_community_matrix.R` | Build site × species abundance matrix | `output/intermediate/community_matrix.rds` |
| `02b_generate_null_communities.R` | Generate 999 Curveball null communities | `output/intermediate/null_communities.rds` |
| `02c_calc_PD.R` | Faith's PD + SES-PD | `output/intermediate/sesPD.csv` |
| `02d_calc_FD.R` | Functional richness (TPD) + SES-FD | `output/intermediate/sesFD.csv` |
| `02e_assemble_diversity_table.R` | Combine TD, PD, FD, µD into final table | `data/processed/diversity_metrics_per_basin.csv` |

### `03_brt_models/` — BRT modelling and variable selection

| Script | Purpose | Key output |
|--------|---------|------------|
| `03a_assemble_modelling_data.R` | Join diversity metrics + environmental predictors | `data/processed/dataprepared.rds` |
| `03b_variable_selection.R` | Pearson correlation + VIF filtering; Ext. Figs. E1–E2 | `data/processed/dataprepared_filtered.rds` |
| `03c_brt_optimize.R` | Hyperparameter grid search (10-fold spatial CV) | `output/intermediate/brt_best_params.csv` |
| `03d_brt_bootstrap.R` | BRT fitting + spatial block bootstrap (1000 iter.) | `output/intermediate/brt_final_results.rds` |
| `03e_brt_validate.R` | Spatial CV validation + null model comparison | `output/intermediate/brt_cv_summary.csv` |
| `03f_brt_figures.R` | Figs. 3–5; Extended Figs. E4–E7 | `output/figures/main/Fig3_*.png`, etc. |

### `04_conservation_overlap/` — Conservation prioritization analysis

| Script | Purpose | Key output |
|--------|---------|------------|
| `04a_conservation_overlap.R` | UpSet plots (Fig. 2), overlap maps (Ext. Fig. E3), conservation plan alignment (Fig. 6) | `output/figures/main/Fig2_*.png`, `Fig6_*.png` |

### `05_figures/` — Remaining manuscript figures

| Script | Purpose | Key output |
|--------|---------|------------|
| `05a_fig1_rank_maps.R` | Diversity rank maps (Fig. 1) | `output/figures/main/Fig1_*.png` |
| `05b_figS4_distribution_maps.R` | Spatial distribution + latitudinal gradient (SI Fig. S4) | `output/figures/si/FigS4_*.png` |

### `utils/` — Shared functions

| File | Content |
|------|---------|
| `curveball.R` | Curveball null model algorithm (Strona et al. 2014) |
| `clean_name.R` | Species name standardisation utility |
| `R_fast_function.R` | Helper functions for TPD calculations |
| `brt_functions.R` | BRT fitting, bootstrap, PDP, Friedman's H helpers |

### `src/`

| File | Content |
|------|---------|
| `fast_tpd_metrics_omp.cpp` | C++ implementation for fast TPD metric computation (OpenMP parallelised) |

---

## Data Sources

Large raster datasets are not bundled with this repository due to file size. Download them before running `01e`–`01g`.

| Dataset | Source | Reference |
|---------|--------|-----------|
| AmazonFish occurrences & subdrainage basins | [DOI: 10.1038/s41597-020-0436-4](https://doi.org/10.1038/s41597-020-0436-4) | Jézéquel et al. (2020) |
| WorldClim v2.1 (Bio, srad, vapr, elevation) | [worldclim.org](https://worldclim.org/data/worldclim.html) | Fick & Hijmans (2017) |
| MODIS NPP (MOD17A3HGF V006) | [NASA LP DAAC](https://lpdaac.usgs.gov/) | Running et al. (2015) |
| Global PET/AET | [CGIAR-CSI](https://www.cgiar-csi.org/data/global-aridity-and-pet-database/) | Trabucco & Zomer (2010) |
| UNH-GRDC Composite Runoff | [compositerunoff.sr.unh.edu](https://www.compositerunoff.sr.unh.edu/) | Fekete et al. (2000) |
| Paleoclimate (CCSM, MIROC, MPI-ESM-P) | Supplementary of Varela et al. (2015) | Varela et al. (2015) |
| River network (SWORD) | Altenau et al. (2021) | Altenau et al. (2021) |
| FishBase taxonomy | [fishbase.org](https://www.fishbase.org) | Froese & Pauly (2025) |
| FISHMORPH trait database | Su et al. (2021); Brosse et al. (2021) | — |

### ArcGIS-derived intermediates (bundled in `data/processed/`)

The following files were generated in ArcGIS Pro and cannot be regenerated from this pipeline alone:

- `slope15.csv` — % basin area with slope ≥ 15° (from WorldClim DEM)
- `soil_types.csv` — Soil type proportions (Harmonized World Soil Database; Caspari & Schad 2022)
- `land_cover.csv` — Land cover proportions (ESA CCI 2017)
- `water_color.csv` — Water color classification (Venticinque et al. 2016)
- `waterfalls.csv` — Waterfall counts (SWORD database; Altenau et al. 2021)
- `pebas_lake.csv` — Pebas Lake extent indicator (Hoorn et al. 2010)
- `jezequel_overlap.csv` — Proactive conservation plan basin lists (Jézéquel et al. 2020)

---

## Dependencies

### R packages

```r
# Core
tidyverse, here, data.table

# Spatial
sf, terra, exactextractr, blockCV, SpatialPack

# Phylogenetics & traits
FishPhyloMaker, rfishbase, PVR, missForest, picante, TPD, ks, adiv

# Modelling
gbm3, car, psych, vegan, paran, Rmisc

# Plotting
ggplot2, ggpubr, cowplot, patchwork, ComplexUpset, VennDiagram,
eulerr, ggcorrplot2, ggplotify, paletteer, scico

# Parallelisation
parallel, Matrix, Rcpp

# Other
iml, ggtext, tictoc, gridExtra
```

### C++ compilation

The file `src/fast_tpd_metrics_omp.cpp` requires a C++ compiler with OpenMP support. It is compiled at runtime via `Rcpp::sourceCpp()`.

## Computational Notes

 - BRT bootstrap (scripts `03d`-`03e`): ~ 5,600 model fits per response variable (1000 influence + 1000 R² + 200 PDP × 14 predictors + 200 Friedman's H). 
 - FD null model (`02d`): 999 permutations with C++ TPD calculations can be accelerated and parallelized via `mclapply` (Unix) or `parLapply` (Windows).
 - Spatial CV: 10-fold block cross-validation using subdrainage basin polygons as spatial blocks (`blockCV::cv_spatial` with `user_blocks`).

## Figure-to-Script Mapping

| Manuscript figure | Generated by |
| --- | --- |
| Fig. 1 (diversity rank maps) | `05a_fig1_rank_maps.R` |
| Fig. 2 (UpSet overlap) | `04a_conservation_overlap.R` |
| Fig. 3 (BRT relative influence) | `03f_brt_figures.R` |
| Fig. 4 (partial dependence plots) | `03f_brt_figures.R` |
| Fig. 5 (Friedman's H) | `03f_brt_figures.R` |
| Fig. 6 (conservation plan overlap) | `04a_conservation_overlap.R` |
| Ext. Fig. E1 (predictor correlation) | `03b_variable_selection.R` |
| Ext. Fig. E2 (response correlation) | `03b_variable_selection.R` |
| Ext. Fig. E3 (spatial overlap maps) | `04a_conservation_overlap.R` |
| Ext. Figs. E4–E7 (2D PDPs) | `03f_brt_figures.R` |
| Ext. Table E1 (BRT parameters) | `03e_brt_validate.R` |
| SI Fig. S1 (taxonomy summary) | `01d_taxonomy_figures.R` |
| SI Fig. S2 (order richness) | `01d_taxonomy_figures.R` |
| SI Fig. S3 (trait completeness) | `01b_impute_traits.R` |
| SI Fig. S4 (distribution + latitudinal gradient) | `05b_figS4_distribution_maps.R` |

## Acknowledgements

Funding: National Natural Science Foundation of China (32470551, 32571868); Chinese Academy of Sciences President's International Fellowship Initiative (2026PD0172); Chinese Academy of Science-Alliance of International Science Organizations (CAS-ANSO) Scholarship for Young Talents.

## Contact

Zi Xun KUA - kuazixun@ihb.ac.cn
Guohuan SU - suguohuan@ihb.ac.cn
