# This script is to analyze flavivirus-host PPI MiST correlations across the following dataset:
# YFV-human (Kenaston et al, this study)
# WNV-human (Li et al, 2019)
# ZIKV-human (Shah et al, 2018)
# DENV-human (Shah et al, 2018)

# Data was compiled by MWK in 2023 from the respective studies
# PSS created separate files for Capsid, NS2B3, NS3, NS4A, and NS5 since these viral baits were consistently studied across all networks
# The goal of this script is to format data so that each prey has an associated MiST score for each virus
# Calculate correlation coefficients for each viral bait
# Bootstrap the correlations to see if it is significant
# Author: PSS, created 20240916, edited 20250107

# Copyright (C) 2025, Shah Lab, UC Davis
# GNU General Public License v3.0
# email: prsshah@ucdavis.edu

# load the dataset
setwd("~/Dropbox/Work/YFV/Conservation/")
library(tidyr)
library(dplyr)

Capsid <- read.csv("Capsid MiST.csv")

# reformat the dataset
Capsid <- Capsid %>% select(Prey,MIST,Virus)
Capsid_wide <- pivot_wider(Capsid, names_from = "Prey", values_from = "MIST")
Capsid_wide <- Capsid_wide %>% select(-Virus)
rownames(Capsid_wide) <- c("DENV", "ZIKVfp", "ZIKVug", "WNV", "YFV")

# create dataframes of pairwise sets and remove any NAs

# DENV vs all
DENVvsZIKVfp <- Capsid_wide[c("DENV", "ZIKVfp"), ]
DENVvsZIKVfp_clean <- DENVvsZIKVfp[, !apply(DENVvsZIKVfp, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVfp_clean) <- c("DENV", "ZIKVfp")
DENVvsZIKVfp_correlation <- cor(as.numeric(DENVvsZIKVfp_clean[1, ]), as.numeric(DENVvsZIKVfp_clean[2, ]), method = "pearson")

DENVvsZIKVug <- Capsid_wide[c("DENV", "ZIKVug"), ]
DENVvsZIKVug_clean <- DENVvsZIKVug[, !apply(DENVvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVug_clean) <- c("DENV", "ZIKVug")
DENVvsZIKVug_correlation <- cor(as.numeric(DENVvsZIKVug_clean[1, ]), as.numeric(DENVvsZIKVug_clean[2, ]), method = "pearson")

DENVvsWNV <- Capsid_wide[c("DENV", "WNV"), ]
DENVvsWNV_clean <- DENVvsWNV[, !apply(DENVvsWNV, 2, function(x) any(is.na(x)))]
rownames(DENVvsWNV_clean) <- c("DENV", "WNV")
DENVvsWNV_correlation <- cor(as.numeric(DENVvsWNV_clean[1, ]), as.numeric(DENVvsWNV_clean[2, ]), method = "pearson")

DENVvsYFV <- Capsid_wide[c("DENV", "YFV"), ]
DENVvsYFV_clean <- DENVvsYFV[, !apply(DENVvsYFV, 2, function(x) any(is.na(x)))]
rownames(DENVvsYFV_clean) <- c("DENV", "YFV")
DENVvsYFV_correlation <- cor(as.numeric(DENVvsYFV_clean[1, ]), as.numeric(DENVvsYFV_clean[2, ]), method = "pearson")

# ZIKVfp vs all
ZIKVfpvsZIKVug <- Capsid_wide[c("ZIKVfp", "ZIKVug"), ]
ZIKVfpvsZIKVug_clean <- ZIKVfpvsZIKVug[, !apply(ZIKVfpvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsZIKVug_clean) <- c("ZIKVfp", "ZIKVug")
ZIKVfpvsZIKVug_correlation <- cor(as.numeric(ZIKVfpvsZIKVug_clean[1, ]), as.numeric(ZIKVfpvsZIKVug_clean[2, ]), method = "pearson")

ZIKVfpvsWNV <- Capsid_wide[c("ZIKVfp", "WNV"), ]
ZIKVfpvsWNV_clean <- ZIKVfpvsWNV[, !apply(ZIKVfpvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsWNV_clean) <- c("ZIKVfp", "WNV")
ZIKVfpvsWNV_correlation <- cor(as.numeric(ZIKVfpvsWNV_clean[1, ]), as.numeric(ZIKVfpvsWNV_clean[2, ]), method = "pearson")

ZIKVfpvsYFV <- Capsid_wide[c("ZIKVfp", "YFV"), ]
ZIKVfpvsYFV_clean <- ZIKVfpvsYFV[, !apply(ZIKVfpvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsYFV_clean) <- c("ZIKVfp", "YFV")
ZIKVfpvsYFV_correlation <- cor(as.numeric(ZIKVfpvsYFV_clean[1, ]), as.numeric(ZIKVfpvsYFV_clean[2, ]), method = "pearson")

# ZIKVug vs all
ZIKVugvsWNV <- Capsid_wide[c("ZIKVug", "WNV"), ]
ZIKVugvsWNV_clean <- ZIKVugvsWNV[, !apply(ZIKVugvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsWNV_clean) <- c("ZIKVug", "WNV")
ZIKVugvsWNV_correlation <- cor(as.numeric(ZIKVugvsWNV_clean[1, ]), as.numeric(ZIKVugvsWNV_clean[2, ]), method = "pearson")

ZIKVugvsYFV <- Capsid_wide[c("ZIKVug", "YFV"), ]
ZIKVugvsYFV_clean <- ZIKVugvsYFV[, !apply(ZIKVugvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsYFV_clean) <- c("ZIKVug", "YFV")
ZIKVugvsYFV_correlation <- cor(as.numeric(ZIKVugvsYFV_clean[1, ]), as.numeric(ZIKVugvsYFV_clean[2, ]), method = "pearson")

# WNV vs all
WNVvsYFV <- Capsid_wide[c("WNV", "YFV"), ]
WNVvsYFV_clean <- WNVvsYFV[, !apply(WNVvsYFV, 2, function(x) any(is.na(x)))]
rownames(WNVvsYFV_clean) <- c("WNV", "YFV")
WNVvsYFV_correlation <- cor(as.numeric(WNVvsYFV_clean[1, ]), as.numeric(WNVvsYFV_clean[2, ]), method = "pearson")

# store all correlations in a table

Capsid_cor <- data.frame(c(1.0, DENVvsZIKVfp_correlation, DENVvsZIKVug_correlation, DENVvsWNV_correlation, DENVvsYFV_correlation), c(DENVvsZIKVfp_correlation, 1.0, ZIKVfpvsZIKVug_correlation, ZIKVfpvsWNV_correlation, ZIKVfpvsYFV_correlation), c(DENVvsZIKVug_correlation, ZIKVfpvsZIKVug_correlation, 1.0, ZIKVugvsWNV_correlation, ZIKVugvsYFV_correlation), c(DENVvsWNV_correlation, ZIKVfpvsWNV_correlation, ZIKVugvsWNV_correlation, 1.0, WNVvsYFV_correlation), c(DENVvsYFV_correlation, ZIKVfpvsYFV_correlation, ZIKVugvsYFV_correlation, WNVvsYFV_correlation, 1.0))

# print dataset size being analyzed
all_objects <- ls()
clean_vectors <- grep("_clean", all_objects, value = TRUE)
for (vec_name in clean_vectors) {
  vec_length <- length(get(vec_name))  # Use get() to retrieve the vector by name
  print(paste("Length of", vec_name, "is", vec_length))
}

# repeat this for NS2B3, NS3, NS4A and NS5

# NS2B3

NS2B3 <- read.csv("NS2B3 MiST.csv")

# reformat the dataset
NS2B3 <- NS2B3 %>% select(Prey,MIST,Virus)
NS2B3_wide <- pivot_wider(NS2B3, names_from = "Prey", values_from = "MIST")
NS2B3_wide <- NS2B3_wide %>% select(-Virus)
rownames(NS2B3_wide) <- c("DENV", "ZIKVfp", "ZIKVug", "WNV", "YFV")

# create dataframes of pairwise sets and remove any NAs

# DENV vs all
DENVvsZIKVfp <- NS2B3_wide[c("DENV", "ZIKVfp"), ]
DENVvsZIKVfp_clean <- DENVvsZIKVfp[, !apply(DENVvsZIKVfp, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVfp_clean) <- c("DENV", "ZIKVfp")
DENVvsZIKVfp_correlation <- cor(as.numeric(DENVvsZIKVfp_clean[1, ]), as.numeric(DENVvsZIKVfp_clean[2, ]), method = "pearson")

DENVvsZIKVug <- NS2B3_wide[c("DENV", "ZIKVug"), ]
DENVvsZIKVug_clean <- DENVvsZIKVug[, !apply(DENVvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVug_clean) <- c("DENV", "ZIKVug")
DENVvsZIKVug_correlation <- cor(as.numeric(DENVvsZIKVug_clean[1, ]), as.numeric(DENVvsZIKVug_clean[2, ]), method = "pearson")

DENVvsWNV <- NS2B3_wide[c("DENV", "WNV"), ]
DENVvsWNV_clean <- DENVvsWNV[, !apply(DENVvsWNV, 2, function(x) any(is.na(x)))]
rownames(DENVvsWNV_clean) <- c("DENV", "WNV")
DENVvsWNV_correlation <- cor(as.numeric(DENVvsWNV_clean[1, ]), as.numeric(DENVvsWNV_clean[2, ]), method = "pearson")

DENVvsYFV <- NS2B3_wide[c("DENV", "YFV"), ]
DENVvsYFV_clean <- DENVvsYFV[, !apply(DENVvsYFV, 2, function(x) any(is.na(x)))]
rownames(DENVvsYFV_clean) <- c("DENV", "YFV")
DENVvsYFV_correlation <- cor(as.numeric(DENVvsYFV_clean[1, ]), as.numeric(DENVvsYFV_clean[2, ]), method = "pearson")

# ZIKVfp vs all
ZIKVfpvsZIKVug <- NS2B3_wide[c("ZIKVfp", "ZIKVug"), ]
ZIKVfpvsZIKVug_clean <- ZIKVfpvsZIKVug[, !apply(ZIKVfpvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsZIKVug_clean) <- c("ZIKVfp", "ZIKVug")
ZIKVfpvsZIKVug_correlation <- cor(as.numeric(ZIKVfpvsZIKVug_clean[1, ]), as.numeric(ZIKVfpvsZIKVug_clean[2, ]), method = "pearson")

ZIKVfpvsWNV <- NS2B3_wide[c("ZIKVfp", "WNV"), ]
ZIKVfpvsWNV_clean <- ZIKVfpvsWNV[, !apply(ZIKVfpvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsWNV_clean) <- c("ZIKVfp", "WNV")
ZIKVfpvsWNV_correlation <- cor(as.numeric(ZIKVfpvsWNV_clean[1, ]), as.numeric(ZIKVfpvsWNV_clean[2, ]), method = "pearson")

ZIKVfpvsYFV <- NS2B3_wide[c("ZIKVfp", "YFV"), ]
ZIKVfpvsYFV_clean <- ZIKVfpvsYFV[, !apply(ZIKVfpvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsYFV_clean) <- c("ZIKVfp", "YFV")
ZIKVfpvsYFV_correlation <- cor(as.numeric(ZIKVfpvsYFV_clean[1, ]), as.numeric(ZIKVfpvsYFV_clean[2, ]), method = "pearson")

# ZIKVug vs all
ZIKVugvsWNV <- NS2B3_wide[c("ZIKVug", "WNV"), ]
ZIKVugvsWNV_clean <- ZIKVugvsWNV[, !apply(ZIKVugvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsWNV_clean) <- c("ZIKVug", "WNV")
ZIKVugvsWNV_correlation <- cor(as.numeric(ZIKVugvsWNV_clean[1, ]), as.numeric(ZIKVugvsWNV_clean[2, ]), method = "pearson")

ZIKVugvsYFV <- NS2B3_wide[c("ZIKVug", "YFV"), ]
ZIKVugvsYFV_clean <- ZIKVugvsYFV[, !apply(ZIKVugvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsYFV_clean) <- c("ZIKVug", "YFV")
ZIKVugvsYFV_correlation <- cor(as.numeric(ZIKVugvsYFV_clean[1, ]), as.numeric(ZIKVugvsYFV_clean[2, ]), method = "pearson")

# WNV vs all
WNVvsYFV <- NS2B3_wide[c("WNV", "YFV"), ]
WNVvsYFV_clean <- WNVvsYFV[, !apply(WNVvsYFV, 2, function(x) any(is.na(x)))]
rownames(WNVvsYFV_clean) <- c("WNV", "YFV")
WNVvsYFV_correlation <- cor(as.numeric(WNVvsYFV_clean[1, ]), as.numeric(WNVvsYFV_clean[2, ]), method = "pearson")

# store all correlations in a table

NS2B3_cor <- data.frame(c(1.0, DENVvsZIKVfp_correlation, DENVvsZIKVug_correlation, DENVvsWNV_correlation, DENVvsYFV_correlation), c(DENVvsZIKVfp_correlation, 1.0, ZIKVfpvsZIKVug_correlation, ZIKVfpvsWNV_correlation, ZIKVfpvsYFV_correlation), c(DENVvsZIKVug_correlation, ZIKVfpvsZIKVug_correlation, 1.0, ZIKVugvsWNV_correlation, ZIKVugvsYFV_correlation), c(DENVvsWNV_correlation, ZIKVfpvsWNV_correlation, ZIKVugvsWNV_correlation, 1.0, WNVvsYFV_correlation), c(DENVvsYFV_correlation, ZIKVfpvsYFV_correlation, ZIKVugvsYFV_correlation, WNVvsYFV_correlation, 1.0))

# print dataset size being analyzed
all_objects <- ls()
clean_vectors <- grep("_clean", all_objects, value = TRUE)
for (vec_name in clean_vectors) {
  vec_length <- length(get(vec_name))  # Use get() to retrieve the vector by name
  print(paste("Length of", vec_name, "is", vec_length))
}

# NS3

NS3 <- read.csv("NS3 MiST.csv")

# reformat the dataset
NS3 <- NS3 %>% select(Prey,MIST,Virus)
NS3_wide <- pivot_wider(NS3, names_from = "Prey", values_from = "MIST")
NS3_wide <- NS3_wide %>% select(-Virus)
rownames(NS3_wide) <- c("DENV", "ZIKVfp", "ZIKVug", "WNV", "YFV")

# create dataframes of pairwise sets and remove any NAs

# DENV vs all
DENVvsZIKVfp <- NS3_wide[c("DENV", "ZIKVfp"), ]
DENVvsZIKVfp_clean <- DENVvsZIKVfp[, !apply(DENVvsZIKVfp, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVfp_clean) <- c("DENV", "ZIKVfp")
DENVvsZIKVfp_correlation <- cor(as.numeric(DENVvsZIKVfp_clean[1, ]), as.numeric(DENVvsZIKVfp_clean[2, ]), method = "pearson")

DENVvsZIKVug <- NS3_wide[c("DENV", "ZIKVug"), ]
DENVvsZIKVug_clean <- DENVvsZIKVug[, !apply(DENVvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVug_clean) <- c("DENV", "ZIKVug")
DENVvsZIKVug_correlation <- cor(as.numeric(DENVvsZIKVug_clean[1, ]), as.numeric(DENVvsZIKVug_clean[2, ]), method = "pearson")

DENVvsWNV <- NS3_wide[c("DENV", "WNV"), ]
DENVvsWNV_clean <- DENVvsWNV[, !apply(DENVvsWNV, 2, function(x) any(is.na(x)))]
rownames(DENVvsWNV_clean) <- c("DENV", "WNV")
DENVvsWNV_correlation <- cor(as.numeric(DENVvsWNV_clean[1, ]), as.numeric(DENVvsWNV_clean[2, ]), method = "pearson")

DENVvsYFV <- NS3_wide[c("DENV", "YFV"), ]
DENVvsYFV_clean <- DENVvsYFV[, !apply(DENVvsYFV, 2, function(x) any(is.na(x)))]
rownames(DENVvsYFV_clean) <- c("DENV", "YFV")
DENVvsYFV_correlation <- cor(as.numeric(DENVvsYFV_clean[1, ]), as.numeric(DENVvsYFV_clean[2, ]), method = "pearson")

# ZIKVfp vs all
ZIKVfpvsZIKVug <- NS3_wide[c("ZIKVfp", "ZIKVug"), ]
ZIKVfpvsZIKVug_clean <- ZIKVfpvsZIKVug[, !apply(ZIKVfpvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsZIKVug_clean) <- c("ZIKVfp", "ZIKVug")
ZIKVfpvsZIKVug_correlation <- cor(as.numeric(ZIKVfpvsZIKVug_clean[1, ]), as.numeric(ZIKVfpvsZIKVug_clean[2, ]), method = "pearson")

ZIKVfpvsWNV <- NS3_wide[c("ZIKVfp", "WNV"), ]
ZIKVfpvsWNV_clean <- ZIKVfpvsWNV[, !apply(ZIKVfpvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsWNV_clean) <- c("ZIKVfp", "WNV")
ZIKVfpvsWNV_correlation <- cor(as.numeric(ZIKVfpvsWNV_clean[1, ]), as.numeric(ZIKVfpvsWNV_clean[2, ]), method = "pearson")

ZIKVfpvsYFV <- NS3_wide[c("ZIKVfp", "YFV"), ]
ZIKVfpvsYFV_clean <- ZIKVfpvsYFV[, !apply(ZIKVfpvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsYFV_clean) <- c("ZIKVfp", "YFV")
ZIKVfpvsYFV_correlation <- cor(as.numeric(ZIKVfpvsYFV_clean[1, ]), as.numeric(ZIKVfpvsYFV_clean[2, ]), method = "pearson")

# ZIKVug vs all
ZIKVugvsWNV <- NS3_wide[c("ZIKVug", "WNV"), ]
ZIKVugvsWNV_clean <- ZIKVugvsWNV[, !apply(ZIKVugvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsWNV_clean) <- c("ZIKVug", "WNV")
ZIKVugvsWNV_correlation <- cor(as.numeric(ZIKVugvsWNV_clean[1, ]), as.numeric(ZIKVugvsWNV_clean[2, ]), method = "pearson")

ZIKVugvsYFV <- NS3_wide[c("ZIKVug", "YFV"), ]
ZIKVugvsYFV_clean <- ZIKVugvsYFV[, !apply(ZIKVugvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsYFV_clean) <- c("ZIKVug", "YFV")
ZIKVugvsYFV_correlation <- cor(as.numeric(ZIKVugvsYFV_clean[1, ]), as.numeric(ZIKVugvsYFV_clean[2, ]), method = "pearson")

# WNV vs all
WNVvsYFV <- NS3_wide[c("WNV", "YFV"), ]
WNVvsYFV_clean <- WNVvsYFV[, !apply(WNVvsYFV, 2, function(x) any(is.na(x)))]
rownames(WNVvsYFV_clean) <- c("WNV", "YFV")
WNVvsYFV_correlation <- cor(as.numeric(WNVvsYFV_clean[1, ]), as.numeric(WNVvsYFV_clean[2, ]), method = "pearson")

# store all correlations in a table

NS3_cor <- data.frame(c(1.0, DENVvsZIKVfp_correlation, DENVvsZIKVug_correlation, DENVvsWNV_correlation, DENVvsYFV_correlation), c(DENVvsZIKVfp_correlation, 1.0, ZIKVfpvsZIKVug_correlation, ZIKVfpvsWNV_correlation, ZIKVfpvsYFV_correlation), c(DENVvsZIKVug_correlation, ZIKVfpvsZIKVug_correlation, 1.0, ZIKVugvsWNV_correlation, ZIKVugvsYFV_correlation), c(DENVvsWNV_correlation, ZIKVfpvsWNV_correlation, ZIKVugvsWNV_correlation, 1.0, WNVvsYFV_correlation), c(DENVvsYFV_correlation, ZIKVfpvsYFV_correlation, ZIKVugvsYFV_correlation, WNVvsYFV_correlation, 1.0))

# print dataset size being analyzed
all_objects <- ls()
clean_vectors <- grep("_clean", all_objects, value = TRUE)
for (vec_name in clean_vectors) {
  vec_length <- length(get(vec_name))  # Use get() to retrieve the vector by name
  print(paste("Length of", vec_name, "is", vec_length))
}

# NS4A

NS4A <- read.csv("NS4A MiST.csv")

# reformat the dataset
NS4A <- NS4A %>% select(Prey,MIST,Virus)
NS4A_wide <- pivot_wider(NS4A, names_from = "Prey", values_from = "MIST")
NS4A_wide <- NS4A_wide %>% select(-Virus)
rownames(NS4A_wide) <- c("DENV", "ZIKVfp", "ZIKVug", "WNV", "YFV")

# create dataframes of pairwise sets and remove any NAs

# DENV vs all
DENVvsZIKVfp <- NS4A_wide[c("DENV", "ZIKVfp"), ]
DENVvsZIKVfp_clean <- DENVvsZIKVfp[, !apply(DENVvsZIKVfp, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVfp_clean) <- c("DENV", "ZIKVfp")
DENVvsZIKVfp_correlation <- cor(as.numeric(DENVvsZIKVfp_clean[1, ]), as.numeric(DENVvsZIKVfp_clean[2, ]), method = "pearson")

DENVvsZIKVug <- NS4A_wide[c("DENV", "ZIKVug"), ]
DENVvsZIKVug_clean <- DENVvsZIKVug[, !apply(DENVvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVug_clean) <- c("DENV", "ZIKVug")
DENVvsZIKVug_correlation <- cor(as.numeric(DENVvsZIKVug_clean[1, ]), as.numeric(DENVvsZIKVug_clean[2, ]), method = "pearson")

DENVvsWNV <- NS4A_wide[c("DENV", "WNV"), ]
DENVvsWNV_clean <- DENVvsWNV[, !apply(DENVvsWNV, 2, function(x) any(is.na(x)))]
rownames(DENVvsWNV_clean) <- c("DENV", "WNV")
DENVvsWNV_correlation <- cor(as.numeric(DENVvsWNV_clean[1, ]), as.numeric(DENVvsWNV_clean[2, ]), method = "pearson")

DENVvsYFV <- NS4A_wide[c("DENV", "YFV"), ]
DENVvsYFV_clean <- DENVvsYFV[, !apply(DENVvsYFV, 2, function(x) any(is.na(x)))]
rownames(DENVvsYFV_clean) <- c("DENV", "YFV")
DENVvsYFV_correlation <- cor(as.numeric(DENVvsYFV_clean[1, ]), as.numeric(DENVvsYFV_clean[2, ]), method = "pearson")

# ZIKVfp vs all
ZIKVfpvsZIKVug <- NS4A_wide[c("ZIKVfp", "ZIKVug"), ]
ZIKVfpvsZIKVug_clean <- ZIKVfpvsZIKVug[, !apply(ZIKVfpvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsZIKVug_clean) <- c("ZIKVfp", "ZIKVug")
ZIKVfpvsZIKVug_correlation <- cor(as.numeric(ZIKVfpvsZIKVug_clean[1, ]), as.numeric(ZIKVfpvsZIKVug_clean[2, ]), method = "pearson")

ZIKVfpvsWNV <- NS4A_wide[c("ZIKVfp", "WNV"), ]
ZIKVfpvsWNV_clean <- ZIKVfpvsWNV[, !apply(ZIKVfpvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsWNV_clean) <- c("ZIKVfp", "WNV")
ZIKVfpvsWNV_correlation <- cor(as.numeric(ZIKVfpvsWNV_clean[1, ]), as.numeric(ZIKVfpvsWNV_clean[2, ]), method = "pearson")

ZIKVfpvsYFV <- NS4A_wide[c("ZIKVfp", "YFV"), ]
ZIKVfpvsYFV_clean <- ZIKVfpvsYFV[, !apply(ZIKVfpvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsYFV_clean) <- c("ZIKVfp", "YFV")
ZIKVfpvsYFV_correlation <- cor(as.numeric(ZIKVfpvsYFV_clean[1, ]), as.numeric(ZIKVfpvsYFV_clean[2, ]), method = "pearson")

# ZIKVug vs all
ZIKVugvsWNV <- NS4A_wide[c("ZIKVug", "WNV"), ]
ZIKVugvsWNV_clean <- ZIKVugvsWNV[, !apply(ZIKVugvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsWNV_clean) <- c("ZIKVug", "WNV")
ZIKVugvsWNV_correlation <- cor(as.numeric(ZIKVugvsWNV_clean[1, ]), as.numeric(ZIKVugvsWNV_clean[2, ]), method = "pearson")

ZIKVugvsYFV <- NS4A_wide[c("ZIKVug", "YFV"), ]
ZIKVugvsYFV_clean <- ZIKVugvsYFV[, !apply(ZIKVugvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsYFV_clean) <- c("ZIKVug", "YFV")
ZIKVugvsYFV_correlation <- cor(as.numeric(ZIKVugvsYFV_clean[1, ]), as.numeric(ZIKVugvsYFV_clean[2, ]), method = "pearson")

# WNV vs all
WNVvsYFV <- NS4A_wide[c("WNV", "YFV"), ]
WNVvsYFV_clean <- WNVvsYFV[, !apply(WNVvsYFV, 2, function(x) any(is.na(x)))]
rownames(WNVvsYFV_clean) <- c("WNV", "YFV")
WNVvsYFV_correlation <- cor(as.numeric(WNVvsYFV_clean[1, ]), as.numeric(WNVvsYFV_clean[2, ]), method = "pearson")

# store all correlations in a table

NS4A_cor <- data.frame(c(1.0, DENVvsZIKVfp_correlation, DENVvsZIKVug_correlation, DENVvsWNV_correlation, DENVvsYFV_correlation), c(DENVvsZIKVfp_correlation, 1.0, ZIKVfpvsZIKVug_correlation, ZIKVfpvsWNV_correlation, ZIKVfpvsYFV_correlation), c(DENVvsZIKVug_correlation, ZIKVfpvsZIKVug_correlation, 1.0, ZIKVugvsWNV_correlation, ZIKVugvsYFV_correlation), c(DENVvsWNV_correlation, ZIKVfpvsWNV_correlation, ZIKVugvsWNV_correlation, 1.0, WNVvsYFV_correlation), c(DENVvsYFV_correlation, ZIKVfpvsYFV_correlation, ZIKVugvsYFV_correlation, WNVvsYFV_correlation, 1.0))

# print dataset size being analyzed
all_objects <- ls()
clean_vectors <- grep("_clean", all_objects, value = TRUE)
for (vec_name in clean_vectors) {
  vec_length <- length(get(vec_name))  # Use get() to retrieve the vector by name
  print(paste("Length of", vec_name, "is", vec_length))
}

# NS5

NS5 <- read.csv("NS5 MiST.csv")

# reformat the dataset
NS5 <- NS5 %>% select(Prey,MIST,Virus)
NS5_wide <- pivot_wider(NS5, names_from = "Prey", values_from = "MIST")
NS5_wide <- NS5_wide %>% select(-Virus)
rownames(NS5_wide) <- c("DENV", "ZIKVfp", "ZIKVug", "WNV", "YFV")

# create dataframes of pairwise sets and remove any NAs

# DENV vs all
DENVvsZIKVfp <- NS5_wide[c("DENV", "ZIKVfp"), ]
DENVvsZIKVfp_clean <- DENVvsZIKVfp[, !apply(DENVvsZIKVfp, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVfp_clean) <- c("DENV", "ZIKVfp")
DENVvsZIKVfp_correlation <- cor(as.numeric(DENVvsZIKVfp_clean[1, ]), as.numeric(DENVvsZIKVfp_clean[2, ]), method = "pearson")

DENVvsZIKVug <- NS5_wide[c("DENV", "ZIKVug"), ]
DENVvsZIKVug_clean <- DENVvsZIKVug[, !apply(DENVvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(DENVvsZIKVug_clean) <- c("DENV", "ZIKVug")
DENVvsZIKVug_correlation <- cor(as.numeric(DENVvsZIKVug_clean[1, ]), as.numeric(DENVvsZIKVug_clean[2, ]), method = "pearson")

DENVvsWNV <- NS5_wide[c("DENV", "WNV"), ]
DENVvsWNV_clean <- DENVvsWNV[, !apply(DENVvsWNV, 2, function(x) any(is.na(x)))]
rownames(DENVvsWNV_clean) <- c("DENV", "WNV")
DENVvsWNV_correlation <- cor(as.numeric(DENVvsWNV_clean[1, ]), as.numeric(DENVvsWNV_clean[2, ]), method = "pearson")

DENVvsYFV <- NS5_wide[c("DENV", "YFV"), ]
DENVvsYFV_clean <- DENVvsYFV[, !apply(DENVvsYFV, 2, function(x) any(is.na(x)))]
rownames(DENVvsYFV_clean) <- c("DENV", "YFV")
DENVvsYFV_correlation <- cor(as.numeric(DENVvsYFV_clean[1, ]), as.numeric(DENVvsYFV_clean[2, ]), method = "pearson")

# ZIKVfp vs all
ZIKVfpvsZIKVug <- NS5_wide[c("ZIKVfp", "ZIKVug"), ]
ZIKVfpvsZIKVug_clean <- ZIKVfpvsZIKVug[, !apply(ZIKVfpvsZIKVug, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsZIKVug_clean) <- c("ZIKVfp", "ZIKVug")
ZIKVfpvsZIKVug_correlation <- cor(as.numeric(ZIKVfpvsZIKVug_clean[1, ]), as.numeric(ZIKVfpvsZIKVug_clean[2, ]), method = "pearson")

ZIKVfpvsWNV <- NS5_wide[c("ZIKVfp", "WNV"), ]
ZIKVfpvsWNV_clean <- ZIKVfpvsWNV[, !apply(ZIKVfpvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsWNV_clean) <- c("ZIKVfp", "WNV")
ZIKVfpvsWNV_correlation <- cor(as.numeric(ZIKVfpvsWNV_clean[1, ]), as.numeric(ZIKVfpvsWNV_clean[2, ]), method = "pearson")

ZIKVfpvsYFV <- NS5_wide[c("ZIKVfp", "YFV"), ]
ZIKVfpvsYFV_clean <- ZIKVfpvsYFV[, !apply(ZIKVfpvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVfpvsYFV_clean) <- c("ZIKVfp", "YFV")
ZIKVfpvsYFV_correlation <- cor(as.numeric(ZIKVfpvsYFV_clean[1, ]), as.numeric(ZIKVfpvsYFV_clean[2, ]), method = "pearson")

# ZIKVug vs all
ZIKVugvsWNV <- NS5_wide[c("ZIKVug", "WNV"), ]
ZIKVugvsWNV_clean <- ZIKVugvsWNV[, !apply(ZIKVugvsWNV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsWNV_clean) <- c("ZIKVug", "WNV")
ZIKVugvsWNV_correlation <- cor(as.numeric(ZIKVugvsWNV_clean[1, ]), as.numeric(ZIKVugvsWNV_clean[2, ]), method = "pearson")

ZIKVugvsYFV <- NS5_wide[c("ZIKVug", "YFV"), ]
ZIKVugvsYFV_clean <- ZIKVugvsYFV[, !apply(ZIKVugvsYFV, 2, function(x) any(is.na(x)))]
rownames(ZIKVugvsYFV_clean) <- c("ZIKVug", "YFV")
ZIKVugvsYFV_correlation <- cor(as.numeric(ZIKVugvsYFV_clean[1, ]), as.numeric(ZIKVugvsYFV_clean[2, ]), method = "pearson")

# WNV vs all
WNVvsYFV <- NS5_wide[c("WNV", "YFV"), ]
WNVvsYFV_clean <- WNVvsYFV[, !apply(WNVvsYFV, 2, function(x) any(is.na(x)))]
rownames(WNVvsYFV_clean) <- c("WNV", "YFV")
WNVvsYFV_correlation <- cor(as.numeric(WNVvsYFV_clean[1, ]), as.numeric(WNVvsYFV_clean[2, ]), method = "pearson")

# store all correlations in a table

NS5_cor <- data.frame(c(1.0, DENVvsZIKVfp_correlation, DENVvsZIKVug_correlation, DENVvsWNV_correlation, DENVvsYFV_correlation), c(DENVvsZIKVfp_correlation, 1.0, ZIKVfpvsZIKVug_correlation, ZIKVfpvsWNV_correlation, ZIKVfpvsYFV_correlation), c(DENVvsZIKVug_correlation, ZIKVfpvsZIKVug_correlation, 1.0, ZIKVugvsWNV_correlation, ZIKVugvsYFV_correlation), c(DENVvsWNV_correlation, ZIKVfpvsWNV_correlation, ZIKVugvsWNV_correlation, 1.0, WNVvsYFV_correlation), c(DENVvsYFV_correlation, ZIKVfpvsYFV_correlation, ZIKVugvsYFV_correlation, WNVvsYFV_correlation, 1.0))

# print dataset size being analyzed
all_objects <- ls()
clean_vectors <- grep("_clean", all_objects, value = TRUE)
for (vec_name in clean_vectors) {
  vec_length <- length(get(vec_name))  # Use get() to retrieve the vector by name
  print(paste("Length of", vec_name, "is", vec_length))
}

# remove ZIKVug from analysis to avoid double-weighting ZIKV dataset

Capsid_cor <- Capsid_cor[-3, -3]
NS2B3_cor <- NS2B3_cor[-3, -3]
NS3_cor <- NS3_cor[-3, -3]
NS4A_cor <- NS4A_cor[-3, -3]
NS5_cor <- NS5_cor[-3, -3]

rownames(Capsid_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
colnames(Capsid_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
rownames(NS2B3_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
colnames(NS2B3_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
rownames(NS3_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
colnames(NS3_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
rownames(NS4A_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
colnames(NS4A_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
rownames(NS5_cor) <- c("DENV", "ZIKV", "WNV", "YFV")
colnames(NS5_cor) <- c("DENV", "ZIKV", "WNV", "YFV")

write.csv(Capsid_cor, "Capsid_MiST_correlation.csv")
write.csv(NS2B3_cor, "NS2B3_MiST_correlation.csv")
write.csv(NS3_cor, "NS3_MiST_correlation.csv")
write.csv(NS4A_cor, "NS4A_MiST_correlation.csv")
write.csv(NS5_cor, "NS5_MiST_correlation.csv")

# Now bootstrap the system to test for significance
# Initialize the number of iterations and a vector to store all the correlations
n_iterations <- 10000
Correlation_distribution <- rep(0,n_iterations)
mean_MiST_correlations <- c(0.420, 0.267, 0.075, 0.247, 0.217)
MiST_correlations <- c(0.391032109,	0.581768049,	0.392377484, 0.553373567,	0.215902796, 0.383195056, 
                       0.157173328,	0.122906147,	0.229721064, 0.105340779,	-0.184741949, 0.01804737,
                       0.31742717,	0.321589088,	0.370519845, 0.176180025,	0.267959585, 0.151186309,
                       0.511272229,	0.261221741,	0.209154034, 0.11590231,	0.254872038,	0.128600956,
                       0.320839767,	0.261216522,	0.171943535, 0.346357762,	0.152265007, 0.051812886)*100
# Create a vector of bait sequence identity 
mean_sequence_identity <- c(31.34, 58.52, 54.47, 38.56, 64.17)
sequence_identity <- c(39.39,	37.37,	18.75, 46.15,	24.74, 21.65,
                       61.85,	57.83, 47.99, 63.54,	47.92, 47.66,
                       66.4,	63.59,	51.22, 66.99,	51.22, 51.7,
                       48.82,	42.06,	31.75, 42.86,	34.92, 30.95,
                       66.85,	67.15,	59.35, 70.07,	60.98, 60.62)

Capsid_sequence_identity <- c(39.39,	37.37,	18.75, 46.15,	24.74, 21.65)
NS2B3_sequence_identity <- c(61.85,	57.83, 47.99, 63.54,	47.92, 47.66)
NS3_sequence_identity <- c(66.4,	63.59,	51.22, 66.99,	51.22, 51.7)
NS4A_sequence_identity <- c(48.82,	42.06,	31.75, 42.86,	34.92, 30.95)
NS5_sequence_identity <- c(66.85,	67.15,	59.35, 70.07,	60.98, 60.62)

Capsid_MiST_correlations <- c(0.391032109,	0.581768049,	0.392377484, 0.553373567,	0.215902796, 0.383195056)*100
NS2B3_MiST_correlations <- c(0.157173328,	0.122906147,	0.229721064, 0.105340779,	-0.184741949, 0.01804737)*100
NS3_MiST_correlations <- c(0.31742717,	0.321589088,	0.370519845, 0.176180025,	0.267959585, 0.151186309)*100
NS4A_MiST_correlations <- c(0.511272229,	0.261221741,	0.209154034, 0.11590231,	0.254872038,	0.128600956)*100
NS5_MiST_correlations <- c(0.320839767,	0.261216522,	0.171943535, 0.346357762,	0.152265007, 0.051812886)*100

Capsid_pearsons <- cor(Capsid_MiST_correlations, Capsid_sequence_identity, method = "pearson")
NS2B3_pearsons <- cor(NS2B3_MiST_correlations, NS2B3_sequence_identity, method = "pearson")
NS3_pearsons <- cor(NS3_MiST_correlations, NS3_sequence_identity, method = "pearson")
NS4A_pearsons <- cor(NS4A_MiST_correlations, NS4A_sequence_identity, method = "pearson")
NS5_pearsons <- cor(NS5_MiST_correlations, NS5_sequence_identity, method = "pearson")

print(Capsid_pearsons)
print(NS2B3_pearsons)
print(NS3_pearsons)
print(NS4A_pearsons)
print(NS5_pearsons)

mean_pearsons <- cor(mean_MiST_correlations, mean_sequence_identity, method = "pearson")
MiST_pearsons <- cor(MiST_correlations, sequence_identity, method = "pearson")

for (i in 1:n_iterations) {
  Correlation_distribution[i] <- cor(sample(mean_MiST_correlations), mean_sequence_identity, method = "pearson")
}

threshold <- mean_pearsons
p <- sum(Correlation_distribution < threshold)/n_iterations
print(p)

for (i in 1:n_iterations) {
  Correlation_distribution[i] <- cor(sample(MiST_correlations), sequence_identity, method = "pearson")
}

threshold <- MiST_pearsons
p <- sum(Correlation_distribution < threshold)/n_iterations
print(p)

# test if mean of Capsid MiST correlations is different from the mean of NS5 MiST correlations
t.test(Capsid_MiST_correlations, NS5_MiST_correlations, alternative = "two.sided", var.equal = FALSE)
# test if distribution of Capsid MiST correlations is different from the distribution of NS5 MiST correlations
correlation_data <- data.frame(
  values = c(Capsid_MiST_correlations, NS4A_MiST_correlations, NS2B3_MiST_correlations, NS3_MiST_correlations, NS5_MiST_correlations),
  group = factor(rep(c("Group1", "Group2", "Group3", "Group4", "Group5"), each = 6))
)
kruskal.test(values~group, data = correlation_data)
pairwise.wilcox.test(correlation_data$values, correlation_data$group, p.adjust.method = "none")
wilcox.test(Capsid_MiST_correlations, NS5_MiST_correlations, alternative = "two.sided")

# Randomly shuffle the MiST scores for one bait
# calculate the correlation coefficient
# Do this for n_iterations

for (i in 1:n_iterations) {
  Correlation_distribution[i] <- cor(sample(mean_MiST_correlations), mean_sequence_identity, method = "pearson")
}

threshold <- -0.62
p <- sum(Correlation_distribution < threshold)/n_iterations
print(p)

# Now bootstrap the system to test for significance for the correlation between HC-PPI data and bait Sequence Identity/RMSD
# Initialize the number of iterations and a vector to store all the correlations
n_iterations <- 100000

# The number of prey analyzed for each bait
Capsid_N <- 74
NS3_N <- 24
NS2B3_N <- 20 
NS4A_N <- 55
NS5_N <- 36

# Create a vector of bait sequence identity and a place to store correlations
Sequence_identity <- c(31.34, 58.52, 54.47, 38.56, 64.17)
RMSD <- c(2.765, 1.623, 2.07, 1.69, 1.0)
Correlation_distribution <- rep(0,n_iterations)

# Randomly distribute N prey across 4 bins for each bait
# Do this for n_iterations

for (i in 1:n_iterations) {
  # Create a vector of zeros for each bait
  Capsid.bootstrap <- rep(0, 4)
  NS3.bootstrap <- rep(0, 4)
  NS2B3.bootstrap <- rep(0, 4)
  NS4A.bootstrap <- rep(0, 4)
  NS5.bootstrap <- rep(0, 4)
  # Distribute the prey for this iteration
  for (j in 1:Capsid_N) {
    index <- sample(1:4, 1)
    Capsid.bootstrap[index] <- Capsid.bootstrap[index] + 1
  }
  for (j in 1:NS3_N) {
    index <- sample(1:4, 1)
    NS3.bootstrap[index] <- NS3.bootstrap[index] + 1
  }
  for (j in 1:NS2B3_N) {
    index <- sample(1:4, 1)
    NS2B3.bootstrap[index] <- NS2B3.bootstrap[index] + 1
  }
  for (j in 1:NS4A_N) {
    index <- sample(1:4, 1)
    NS4A.bootstrap[index] <- NS4A.bootstrap[index] + 1
  }
  for (j in 1:NS5_N) {
    index <- sample(1:4, 1)
    NS5.bootstrap[index] <- NS5.bootstrap[index] + 1
  }
  # Calculate the normalized PPI conservation for each bait
  Capsid.avg.PPI <- (Capsid.bootstrap[1]*4+Capsid.bootstrap[2]*3+Capsid.bootstrap[3]*2+Capsid.bootstrap[4]*1)/Capsid_N/4*100
  NS3.avg.PPI <- (NS3.bootstrap[1]*4+NS3.bootstrap[2]*3+NS3.bootstrap[3]*2+NS3.bootstrap[4]*1)/NS3_N/4*100
  NS2B3.avg.PPI <- (NS2B3.bootstrap[1]*4+NS2B3.bootstrap[2]*3+NS2B3.bootstrap[3]*2+NS2B3.bootstrap[4]*1)/NS2B3_N/4*100
  NS4A.avg.PPI <- (NS4A.bootstrap[1]*4+NS4A.bootstrap[2]*3+NS4A.bootstrap[3]*2+NS4A.bootstrap[4]*1)/NS4A_N/4*100
  NS5.avg.PPI <- (NS5.bootstrap[1]*4+NS5.bootstrap[2]*3+NS5.bootstrap[3]*2+NS5.bootstrap[4]*1)/NS5_N/4*100
  PPI_conservation <- c(Capsid.avg.PPI, NS3.avg.PPI, NS2B3.avg.PPI, NS4A.avg.PPI, NS5.avg.PPI)
  # Calculate the correlation between sequence identity of viral baits and conservation of PPIs
  bootstrap.correlation <- cor(Sequence_identity, PPI_conservation, method = "pearson")
  # Store this correlation in a vector to create a distribution of correlations
  Correlation_distribution[i] <- bootstrap.correlation
}

threshold <- -0.89
p <- sum(Correlation_distribution < threshold)/n_iterations
print(p)

for (i in 1:n_iterations) {
  # Create a vector of zeros for each bait
  Capsid.bootstrap <- rep(0, 4)
  NS3.bootstrap <- rep(0, 4)
  NS2B3.bootstrap <- rep(0, 4)
  NS4A.bootstrap <- rep(0, 4)
  NS5.bootstrap <- rep(0, 4)
  # Distribute the prey for this iteration
  for (j in 1:Capsid_N) {
    index <- sample(1:4, 1)
    Capsid.bootstrap[index] <- Capsid.bootstrap[index] + 1
  }
  for (j in 1:NS3_N) {
    index <- sample(1:4, 1)
    NS3.bootstrap[index] <- NS3.bootstrap[index] + 1
  }
  for (j in 1:NS2B3_N) {
    index <- sample(1:4, 1)
    NS2B3.bootstrap[index] <- NS2B3.bootstrap[index] + 1
  }
  for (j in 1:NS4A_N) {
    index <- sample(1:4, 1)
    NS4A.bootstrap[index] <- NS4A.bootstrap[index] + 1
  }
  for (j in 1:NS5_N) {
    index <- sample(1:4, 1)
    NS5.bootstrap[index] <- NS5.bootstrap[index] + 1
  }
  # Calculate the normalized PPI conservation for each bait
  Capsid.avg.PPI <- (Capsid.bootstrap[1]*4+Capsid.bootstrap[2]*3+Capsid.bootstrap[3]*2+Capsid.bootstrap[4]*1)/Capsid_N/4*100
  NS3.avg.PPI <- (NS3.bootstrap[1]*4+NS3.bootstrap[2]*3+NS3.bootstrap[3]*2+NS3.bootstrap[4]*1)/NS3_N/4*100
  NS2B3.avg.PPI <- (NS2B3.bootstrap[1]*4+NS2B3.bootstrap[2]*3+NS2B3.bootstrap[3]*2+NS2B3.bootstrap[4]*1)/NS2B3_N/4*100
  NS4A.avg.PPI <- (NS4A.bootstrap[1]*4+NS4A.bootstrap[2]*3+NS4A.bootstrap[3]*2+NS4A.bootstrap[4]*1)/NS4A_N/4*100
  NS5.avg.PPI <- (NS5.bootstrap[1]*4+NS5.bootstrap[2]*3+NS5.bootstrap[3]*2+NS5.bootstrap[4]*1)/NS5_N/4*100
  PPI_conservation <- c(Capsid.avg.PPI, NS3.avg.PPI, NS2B3.avg.PPI, NS4A.avg.PPI, NS5.avg.PPI)
  # Calculate the correlation between sequence identity of viral baits and conservation of PPIs
  bootstrap.correlation <- cor(RMSD, PPI_conservation, method = "pearson")
  # Store this correlation in a vector to create a distribution of correlations
  Correlation_distribution[i] <- bootstrap.correlation
}

# Plot the histogram of correlations

library(ggplot2)
Correlation_data <- data.frame(Value = Correlation_distribution)
ggplot(Correlation_data, aes(x = Value)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(title = "Distribution of correlation values", x = "Pearson's Correlation", y = "Frequency") +
  theme_minimal()

threshold <- 0.75
p <- sum(Correlation_distribution > threshold)/n_iterations
print(p)
