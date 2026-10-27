# =============================================================================
# Script:  01a_make_phylo_tree.R
# Project: Amazon multidimensional biodiversity
# Purpose: Construct a time-calibrated phylogenetic tree for 2,492 Amazonian
#          freshwater fish species using FishPhyloMaker.
#
# Input:   data/processed/trait_matrix_2492.csv
# Output:  output/derived/phylogeny_amazon_2492.rds
#
# Dependencies: FishPhyloMaker, rfishbase, tidyverse, here
# Notes:   Freshwater stingrays (Potamotrygonidae) are excluded because they
#          are not represented in the FishPhyloMaker backbone.
# =============================================================================
# =============================================================================
# 01a_make_phylo_tree.R
# =============================================================================

# load libraries
library(FishPhyloMaker)
library(tidyverse)
library(here)

# set paths
input_dir <- here::here("data", "processed")
inter_dir <- here::here("output", "intermediate")
output_dir <- here::here("output", "derived")

dir.create(inter_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# load trait data
trait_raw <- read.csv(file.path(input_dir,"trait_miss_2492.csv"))

# load and cache FishBase taxonomy
tax_cache <- here::here("data", "processed", "fishbase_taxonomy_cache.csv")

if (file.exists(tax_cache)) {
  fb_taxa <- read.csv(tax_cache)
} else {
  fb_taxa_raw <- rfishbase::load_taxa(server = "fishbase")
  fb_taxa <- fb_taxa_raw |>
    select(Class, Order, Family, Genus) |>
    distinct()
  write.csv(fb_taxa, tax_cache, row.names = FALSE)
}

# merge taxonomy with trait data
trait <- trait_raw |>
  mutate(Genus = word(Genus.species, 1),
         Species = word(Genus.species, 2)) |>
  left_join(fb_taxa, by = 'Genus') |>
  # fill in two remaining entries for unresolved species
  mutate(Family = case_when(Genus %in% c("Araichthys", "Loraxichthys","Lipopterichthys","Aposturisoma", "Etsaputu") ~ "Loricariidae",
                            Genus == "Melanosternarchus" ~ "Apteronotidae",
                            Genus == "Stegostenopos" ~ "Rhamphichthyidae",
                            Genus == "Piabucina" ~ "Lebiasinidae",
                            TRUE ~ Family),
         Order = case_when(Family == "Loricariidae" ~ "Siluriformes",
                           Family == "Apteronotidae" ~ "Gymnotiformes",
                           Family == "Rhamphichthyidae" ~ "Gymnotiformes",
                           Family == "Lebiasinidae" ~ "Characiformes",
                           TRUE ~ Order),
         Class = case_when(Order == "Siluriformes" ~ "Teleostei",
                           Order == "Gymnotiformes" ~ "Teleostei",
                           Order == "Characiformes" ~ "Teleostei",
                           TRUE ~ Class),
         MBl = ifelse(Genus.species == "Melanosternarchus amaru", 27.1, MBl),
         SpecCode = ifelse(Genus.species == "Melanosternarchus amaru", 71060, SpecCode))

# verify no missing taxonomy
stopifnot(!any(is.na(trait$Family) | is.na(trait$Order) | is.na(trait$Class)))

# build taxon data (uses interactive insertions)
taxon_data <- FishTaxaMaker(str_replace_all(trait$Genus.species," ","_"), allow.manual.insert = TRUE)
Doradidae
Siluriformes
Trichomycteridae
Siluriformes
Rivulidae
Cyprinodontiformes
Rivulidae
Cyprinodontiformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Callichthyidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Iguanodectidae
Characiformes
Iguanodectidae
Characiformes
Iguanodectidae
Characiformes
Iguanodectidae
Characiformes
Heptapteridae
Siluriformes
Heptapteridae
Siluriformes
Heptapteridae
Siluriformes
Crenuchidae
Characiformes
Crenuchidae
Characiformes
Crenuchidae
Characiformes
Stevardiidae
Characiformes
Callichthyidae
Siluriformes
Callichthyidae
Siluriformes
Callichthyidae
Siluriformes
Stevardiidae
Characiformes
Stevardiidae
Characiformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Curimatidae
Characiformes
Curimatidae
Characiformes
Curimatidae
Characiformes
Curimatidae
Characiformes
Curimatidae
Characiformes
Curimatidae
Characiformes
Sternopygidae
Gymnotiformes
Sternopygidae
Gymnotiformes
Characidae
Characiformes
Gymnotidae
Gymnotiformes
Gymnotidae
Gymnotiformes
Gymnotidae
Gymnotiformes
Gymnotidae
Gymnotiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Hemiodontidae
Characiformes
Hemiodontidae
Characiformes
Erythrinidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Heptapteridae
Siluriformes
Stevardiidae
Characiformes
Stevardiidae
Characiformes
Stevardiidae
Characiformes
Rivulidae
Cyprinodontiformes
Rivulidae
Cyprinodontiformes
Rivulidae
Cyprinodontiformes
Anostomidae
Characiformes
Heptapteridae
Siluriformes
Rivulidae
Cyprinodontiformes
Apteronotidae
Gymnotiformes
Serrasalmidae
Characiformes
Crenuchidae
Characiformes
Eleotridae
Gobiiformes
Characidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Serrasalmidae
Characiformes
Loricariidae
Siluriformes
Rivulidae
Cyprinodontiformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Loricariidae
Siluriformes
Stevardiidae
Characiformes
Acestrorhamphidae
Characiformes
Acestrorhamphidae
Characiformes
Stevardiidae
Characiformes
Stevardiidae
Characiformes
Stevardiidae
Characiformes
Cichlidae
Cichliformes
Anostomidae
Characiformes
Rivulidae
Cyprinodontiformes
Auchenipteridae
Siluriformes
Auchenipteridae
Siluriformes
Acestrorhamphidae
Characiformes
Trichomycteridae
Siluriformes
Trichomycteridae
Siluriformes
Acestrorhamphidae
Characiformes
Iguanodectidae
Characiformes
Crenuchidae
Characiformes
Rivulidae
Cyprinodontiformes
Rivulidae
Cyprinodontiformes
Eleotridae
Gobiiformes

# save taxon_data
saveRDS(taxon_data, paste0(path, "2b.phylo_taxon_data.rds"))

# resolve remaining unmatched species
taxdata1 <- taxon_data$All_info_fishbase |>
  mutate(valid_names = case_when(user_spp == "Apteronotus_quilombola" ~ "Apteronotus quilombola",
                                 user_spp == "Araichthys_loro" ~ "Araichthys loro",
                                 TRUE ~ valid_names),
         Family = case_when(valid_names == "Apteronotus quilombola" ~ "Apteronotidae",
                            valid_names == "Araichthys loro" ~ "Loricariidae",
                            TRUE ~ Family),
         Order = case_when(Family == "Apteronotidae" ~ "Gymnotiformes",
                           Family == "Loricariidae" ~ "Siluriformes",
                           TRUE ~ Order),
         Class = case_when(Order == "Gymnotiformes" ~ "Teleostei",
                           Order == "Siluriformes" ~ "Teleostei",
                           TRUE ~ Class)) |>
  rownames_to_column(var = "SpecCode")
# get tax data to be used for building phylogenetic tree
taxdata2 <- taxon_data$Taxon_data_FishPhyloMaker |>
  mutate(s = ifelse(s == "Genus_Species","Apteronotus_quilombola",s),
         f = ifelse(s == "Apteronotus_quilombola","Apteronotidae",f),
         o = ifelse(o == "Apteronotidae","Gymnotiformes",o)) |>
  drop_na() |> # drop all NAs
  distinct()

# sanity check to make sure no more missing data
stopifnot(sum(is.na(taxdata2$s)) == 0)
stopifnot(sum(is.na(taxdata2$f)) == 0)
stopifnot(sum(is.na(taxdata2$o)) == 0)

# build phylogenetic tree
trait_taxon <- trait |>
  dplyr::select(s = Genus.species, f = Family, o = Order) |>
  mutate(s = str_replace_all(s," ","_"))

## using interactive manual inserts
phylo_2492 <- FishPhyloMaker(trait_taxon, insert.base.node = FALSE, return.insertions = TRUE, progress.bar = FALSE)
Leporinus
Pseudobunocephalus
Myleus
Pimelodus
Aphyocharax
Characidium
Sarcoglanis
Peckoltia
Acanthicus
Anostomus
Hypoclinemus
Farlowella
Corymbophanes
Liosomadoras
Tetragonopterus
Parancistrus
Ctenolucius
Gladioglanis
Pimelodella
Exodon
Copella
Chaetostoma
Corumbataia
Potamorhina
Prionobrama
Cetopsis
Archolaemus
Tatia
Ammocryptocharax
Gymnotus
Ageneiosus
Micromyzon
Peckoltia
Pseudacanthicus
Centromochlus
Loricariichthys
Ferrarissoaresia
Characidium
Vandellia
Anostomus
Araichthys
Loricariichthys
Araichthys
Imparfinis
Leporinus
Hypopomidae
Rhamphichthyidae
Trichomycterus
Elachocharax
Anostomoides
Nannostomus
Anostomidae
Characidium
Eleotridae
Brachyglanis
Hypoptopoma
Loricariichthys
Chaetostoma
Gilbertolus
Chaetostoma
Crenicichla
Sarcoglanis
Imparfinis
Apteronotus
Abramites
Ammocryptocharax
Compsaraia
Elachocharax
Anodus
Vandellia
Gladioglanis
Nemuroglanis
Hypoptopoma
Hypoptopoma
Dekeyseria
Acestridium
Microcharacidium
Plagioscion
Plagioscion
Vandellia
Chaetostoma
Cetopsorhamdia
Compsaraia
Pachyurus
Anostomus
Nannostomus
Pseudohemiodon
Pariosternarchus
Vandellia
Sternarchogiton
Sarcoglanis
Ilisha
Microsternarchus
Exallodontus
Myloplus
Leptagoniates
Petulanos
Asterophysus
Loricaria
Rhinelepis
Sturisoma
Typhlobelus
Brachyhypopomus
Harttia
Apistoloricaria
Goeldiella
Parotocinclus
Pseudopimelodus
Schizodon
Gnathodolus
Homodiaetus
Peckoltia
Leporacanthicus
Simpsonichthys
Ageneiosus
Sarcoglanis
Hypopygus
Ochmacanthus
Hypopygus
Stauroglanis
Sternarchogiton
Platyurosternarchus
Anostomus
Monopterus
Nemadoras
Ageneiosus
Liosomadoras
Auchenipterichthys
Auchenipterichthys
Miuroglanis
Tridens
Tridens
Ageneiosus
Myloplus
Crenicichla
Pterobunocephalus
Pimelodidae
Characidae
Characidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Loricariidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Poeciliidae
Poeciliidae
Poeciliidae
Poeciliidae
Poeciliidae
Poeciliidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characiformes
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characiformes
Characidae
Erythrinidae
Characiformes
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae
Characidae

# save tree
saveRDS(phylo_2492, file.path(output_dir, "phylogeny_amazon_2492.rds"))
message("Phylogenetic tree saved: ", nrow(taxdata2) - nrow(species_miss), " tips")

# identify excluded species
miss_index <- match(setdiff(taxdata2$s, phylo_2492$Phylogeny$tip.label), taxdata2$s)
species_miss <- data.frame(s=taxdata2$s[miss_index], f=taxdata2$f[miss_index], o=taxdata2$o[miss_index])
rownames(species_miss) <- miss_index

message("Species excluded from phylogeny: ", nrow(species_miss))
print(species_miss)

# Final validation
tree <- phylo_2492$Phylogeny
message("Tree tips: ", length(tree$tip.label))
message("Expected species: ", nrow(taxdata2) - nrow(species_miss))
stopifnot(length(tree$tip.label) == nrow(taxdata2) - nrow(species_miss))