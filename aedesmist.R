### Aedes aegypti MiST Scoring and Human Homolog Mapping
### Description: MiST scoring pipeline for Aedes aegypti proteomics data, with mapping to human homologs using InParanoid.
### The MiST package was developed by the Krogan Lab at UCSF: https://doi.org/10.1002/0471250953.bi0819s49.
### Author: Matthew Kenaston, Date: 23 JUN 2025

### Copyright (C) 2025, Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu

#### LOAD REQUIRED LIBRARIES ####
library(MiST.ShahLab)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

#### STEP 1: CLEAN AND FORMAT RAW PROTEIN GROUPS ####
raw <- read.csv("data/aeco_raw.csv")

# Filter for unique protein groups (no semicolon = non-redundant)
unique_rows <- grep(";", raw$PG.UniProtIds, invert = TRUE)
sep_single <- raw[unique_rows, ]

# Expand multi-mapped rows
sep_multi <- raw %>% separate_rows(PG.UniProtIds, sep = ";")
write.csv(sep_multi, "output/aeco_pg_separated.csv")

#### STEP 2: PREPROCESS INPUT FOR MiST ####
input_file <- preprocess.main(
  data_file = "data/input/data_bgrm.txt",
  keys_file = "data/input/keys.txt",
  remove_file = "data/input/remove.txt",
  collapse_file = "data/input/collapse.txt",
  exclusions_file = "data/input/specificity_exclusions.txt",
  output_file = "output/preprocessed_aeco.txt",
  contaminants_file = "data/input/contaminants.txt",
  id_colname = "id",
  prey_colname = "uniprot",
  pepcount_colname = "ms_unique",
  mw_colname = "ms_protein_mw",
  filter_data = 0,
  rm_co = 0
)

# Optional QC step
qc.main(matrix_file = input_file, cluster = 1, font_scale = 9)

#### STEP 3: RUN MiST SCORING ####
mist_result <- mist.main(
  matrix_file = input_file,
  weights = "fixed",
  w_R = 0.35,
  w_S = 0.57,
  w_A = 0.08
)
write.table(mist_result, gsub(".txt", "_MIST.txt", input_file), row.names = FALSE, sep = "	")

#### STEP 4: COMBINE WITH SAINT RESULTS ####
saint_data <- read.csv("data/saint_aeco.csv")
combined <- inner_join(mist_result, saint_data)
write.csv(combined, "output/mist_saint_combined.csv")

# Filter for high-confidence interactions
filtered <- filter(combined, MIST > 0.65 & BFDR <= 0.05)

#### STEP 5: SEPARATE MULTI-GENE ANNOTATIONS AND MAP TO HUMAN ####
split_result <- filtered %>% 
  separate_rows(PreyGene, sep = ";") %>%
  mutate(PreyGene = gsub("_[^_]+$", "", PreyGene))

# Load InParanoid mapping (Aedes-to-Human)
inparanoid_map <- read.csv("data/aedes_to_human_inparanoid.csv") %>%
  filter(score > 0)

# Map each Aedes protein to its next human homolog (highest confidence)
map_to_human <- function(protein) {
  match_index <- which(inparanoid_map$uniprot == protein)
  if (length(match_index) == 0) return("")
  candidate <- inparanoid_map[match_index + 1, ]
  while (nrow(candidate) > 0 && candidate$species == "A.aegypti") {
    match_index <- match_index + 1
    candidate <- inparanoid_map[match_index + 1, ]
  }
  return(ifelse(nrow(candidate) > 0, candidate$uniprot, ""))
}

split_result$h.conv <- sapply(split_result$PreyGene, map_to_human)
split_result$pair <- paste(split_result$Bait, split_result$Prey, sep = "_")

# Keep mapped or fallback to unmapped if unique
mapped <- split_result[split_result$h.conv != "" & !duplicated(split_result$pair), ]
unmapped <- split_result[split_result$h.conv == "" & !duplicated(split_result$pair), ]
final_result <- bind_rows(mapped, unmapped) %>%
  distinct(pair, .keep_all = TRUE)

#### STEP 6: EXPORT MAPPED RESULTS ####
write.csv(split_result, "output/aeco_interolog_comparative_full.csv")
