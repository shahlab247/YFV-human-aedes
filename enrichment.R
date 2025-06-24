### Gene Enrichment Analysis for MiST-Identified Preys
### This script performs GO and pathway enrichment using gProfiler for each viral bait.
### Matthew Kenaston, 23 JUN 2025

### Copyright (C) 2025, Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu

#### LOAD LIBRARIES ####
library(gprofiler2)      # Gene enrichment (GO, KEGG, Reactome)
library(dplyr)           # Data manipulation
library(readr)           # Export to CSV

#### INPUT: CLEANED MiST + SAINT RESULTS ####
# Replace with final filtered results including UniProt IDs
final_mist <- read.csv("data/MIST_final_filtered.csv")

#### PERFORM ENRICHMENT PER BAIT ####
unique_baits <- unique(final_mist$Bait)
enriched_results <- data.frame()

for (bait in unique_baits) {
  bait_preys <- filter(final_mist, Bait == bait)

  enrichment <- gost(
    query = bait_preys$Prey,
    correction_method = "bonferroni",
    domain_scope = "known",
    evcodes = TRUE,
    significant = FALSE
  )$result

  enrichment$Bait <- bait
  enriched_results <- bind_rows(enriched_results, enrichment)
}

#### EXPORT RESULTS ####
write_csv(enriched_results, "output/enrichment_results.csv")
