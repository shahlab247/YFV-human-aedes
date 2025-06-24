### YFV-Human MiST Interaction Analysis Pipeline
### General pipeline for MS-based PPI data scoring using MiST, background filtering, and CORUM complex annotation.
### Matthew Kenaston, 23 JUN 2025

### Copyright (C) 2025, Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu

#### LOAD REQUIRED LIBRARIES ####
library(MiST.ShahLab)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(pheatmap)
library(biomaRt)

#### STEP 1: LOAD & RESHAPE RAW DATA ####
raw <- read.csv("data/input_spectronaut.csv")
base <- rep(raw[, 1:6], 36)
reps <- lapply(1:length(raw[, 7:42]), function(i) {
  data.frame(id = colnames(raw[(6 + i)]), ms_unique = raw[, (6 + i)])
})
reps.all <- do.call("rbind", reps)
predata <- cbind(base, reps.all)
write.csv(predata, "output/predata_unprocessed.csv")

#### STEP 2: PILOT MiST SCORING ####
data <- read.delim("data/input/data.txt") %>%
  filter(ms_unique > 0)
write.table(data, "data/input/data-thresholded.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '	')

input <- preprocess.main(
  data_file = "data/input/data-thresholded.txt",
  keys_file = "data/input/keys.txt",
  remove_file = "data/input/remove.txt",
  collapse_file = "data/input/collapse.txt",
  exclusions_file = "data/input/exclusions.txt",
  output_file = "output/preprocessed_thresholded.txt",
  contaminants_file = "data/input/contaminants.txt",
  id_colname = "id",
  prey_colname = "uniprot",
  pepcount_colname = "ms_unique",
  mw_colname = "ms_protein_mw",
  filter_data = 1,
  rm_co = 0
)

pilot_result <- mist.main(
  matrix_file = input,
  training_file = "data/input/training.txt",
  weights = "fixed",
  w_R = 0.35, w_S = 0.57, w_A = 0.08
)
write.table(pilot_result, "output/MIST_pilot.txt", row.names = FALSE, col.names = TRUE, sep = "	")

#### STEP 3: BACKGROUND REMOVAL BASED ON SAINT ####
saint <- read.csv("data/saint_output.csv")
saint$pair <- paste(saint$Bait, saint$Prey)
pilot_result$pair <- paste(pilot_result$Bait, pilot_result$Prey)
pilot_result$multi <- stringr::str_count(pilot_result$Ip, ".R")
pilot_3hit <- filter(pilot_result, multi > 0)

yfv_baits <- unique(pilot_result$Bait)

# Compute FC relative to prey background mean
fc_table <- data.frame()
for (prey in unique(saint$Prey)) {
  prey_entries <- filter(saint, Prey == prey)
  count <- nrow(prey_entries)
  miss <- length(yfv_baits) - count
  rep_zeros <- rep(0, miss)
  gfp_ctrl <- mean(as.numeric(unlist(strsplit(prey_entries$ctrlCounts[1], split = "|", fixed = TRUE))))
  avg <- mean(c(prey_entries$AvgSpec, rep_zeros, gfp_ctrl))
  prey_entries$FC <- prey_entries$AvgSpec / avg
  fc_table <- rbind(fc_table, prey_entries)
}

# ROC thresholding
library(dismo)
positive <- filter(fc_table, SaintScore >= 0.95 & Prey %in% filter(pilot_result, MIST > 0.67)$Prey)
roc_eval <- evaluate(p = positive$FC, a = fc_table$FC)
fc_thresh <- threshold(roc_eval)$equal_sens_spec

# Filter raw input for final MiST run
data$pair <- paste(gsub("\..*", "", data$id), data$uniprot)
filtered_data <- filter(data, pair %in% filter(fc_table, FC > fc_thresh)$pair & pair %in% pilot_3hit$pair)
write.table(filtered_data, "data/input/data-final.txt", row.names = FALSE, col.names = TRUE, sep = "	", quote = FALSE)

#### STEP 4: FINAL MIST SCORING ####
input_final <- preprocess.main(
  data_file = "data/input/data-final.txt",
  keys_file = "data/input/keys.txt",
  remove_file = "data/input/remove.txt",
  collapse_file = "data/input/collapse.txt",
  exclusions_file = "data/input/exclusions.txt",
  output_file = "output/preprocessed_final.txt",
  contaminants_file = "data/input/contaminants.txt",
  id_colname = "id",
  prey_colname = "uniprot",
  pepcount_colname = "ms_unique",
  mw_colname = "ms_protein_mw",
  filter_data = 0,
  rm_co = 0
)

results_final <- mist.main(
  matrix_file = input_final,
  training_file = "data/input/training.txt",
  weights = "fixed",
  w_R = 0.35, w_S = 0.57, w_A = 0.08
)

write.table(results_final, "output/MIST_final_results.txt", row.names = FALSE, col.names = TRUE, sep = "	")

#### STEP 5: CORUM COMPLEX RESCUE ####
corum <- read.delim("data/corum.txt")
results_final$pair <- paste(results_final$Bait, results_final$Prey)
filtered_final <- filter(results_final, MIST >= 0.6)

rescued_interactions <- data.frame()
for (b in unique(filtered_final$Bait)) {
  bait_set <- filter(filtered_final, Bait == b)
  sig_hits <- filter(bait_set, MIST >= 0.67)
  available_hits <- bait_set$Prey

  corum_hits <- corum[sapply(strsplit(corum$subunits.UniProt.IDs., ";"), function(members) {
    sum(members %in% available_hits) >= 2
  }), ]

  for (i in 1:nrow(corum_hits)) {
    members <- strsplit(corum_hits$subunits.UniProt.IDs.[i], ";")[[1]]
    rescue <- filter(bait_set, Prey %in% members)
    rescued_interactions <- rbind(rescued_interactions, rescue)
  }

  final_bait_set <- rbind(sig_hits, rescued_interactions)
  write.csv(final_bait_set, paste0("output/complex_rescue_", b, ".csv"), row.names = FALSE)
}
