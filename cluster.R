### Flavivirus-Human Interaction Network Analysis with DIS Scores and Enrichment]
### This script performs clustering, calculates virus-specific distinctiveness (DIS) scores, 
### and conducts GO/CORUM/Reactome enrichment for MiST interaction data across flaviviruses.
### Matthew Kenaston, 23 JUN 2025

### Copyright (C) 2025, Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu

#### LOAD REQUIRED LIBRARIES ####
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(cluster)
library(gprofiler2)

#### STEP 1: LOAD FLAVIVIRUS-WIDE MiST DATA ####
mist <- read.csv("data/Flavivirus_Global_MIST.csv")

#### STEP 2: PREPROCESS DATA FOR CLUSTERING ####
mist$pair <- paste(mist$Bait, mist$Prey, sep = "_")
mist$sig <- mist$MIST > 0.67
wide_mist <- pivot_wider(mist, names_from = Virus, values_from = MIST, values_fill = 0)
wide_mist$sig_sum <- rowSums(wide_mist[grep("^sig_", colnames(wide_mist))], na.rm = TRUE)
wide_mist <- wide_mist[wide_mist$sig_sum > 0, ]

#### STEP 3: CLUSTERING ####
mist_matrix <- as.matrix(wide_mist[grep("^MIST_", colnames(wide_mist))])
row.names(mist_matrix) <- wide_mist$pair
hclust_res <- hclust(dist(mist_matrix), method = "ward.D2")
kmeans_res <- kmeans(dist(mist_matrix), centers = 17)

wide_mist$hclust <- cutree(hclust_res, k = 17)
wide_mist$kmeans <- kmeans_res$cluster
wide_mist$Bait <- sub("_[^_]+$", "", wide_mist$pair)
wide_mist$Prey <- sub("^[^_]*_", "", wide_mist$pair)

#### STEP 4: CALCULATE VIRUS-SPECIFIC DIS SCORES ####
# DIS = MiST_V - mean(MiST_others)
wide_mist$MIST_avg <- rowMeans(wide_mist[grep("^MIST_", colnames(wide_mist))])

wide_mist$DIS_YFV  <- wide_mist$MIST_YFV  - rowMeans(wide_mist[, c("MIST_DENV", "MIST_WNV", "MIST_ZIKVug", "MIST_ZIKVfp")])
wide_mist$DIS_DENV <- wide_mist$MIST_DENV - rowMeans(wide_mist[, c("MIST_YFV", "MIST_WNV", "MIST_ZIKVug", "MIST_ZIKVfp")])
wide_mist$DIS_WNV  <- wide_mist$MIST_WNV  - rowMeans(wide_mist[, c("MIST_YFV", "MIST_DENV", "MIST_ZIKVug", "MIST_ZIKVfp")])
wide_mist$DIS_ZIKV <- rowMeans(wide_mist[, c("MIST_ZIKVug", "MIST_ZIKVfp")]) -
                      rowMeans(wide_mist[, c("MIST_YFV", "MIST_DENV", "MIST_WNV")])
wide_mist$DIS_avg <- rowMeans(wide_mist[, c("DIS_YFV", "DIS_DENV", "DIS_WNV", "DIS_ZIKV")])

#### STEP 5: EXPORT DIS SCORES ####
write_csv(wide_mist, "output/flavi_net_clusters_with_DIS.csv")

#### STEP 6: ENRICHMENT ANALYSIS BY CLUSTER AND BAIT ####
wide_mist$group <- paste(wide_mist$Bait, wide_mist$hclust, sep = "_")
enrich_input <- split(wide_mist$Prey, wide_mist$group)

enrich_results <- data.frame()
for (grp in names(enrich_input)) {
  go_res <- gost(enrich_input[[grp]],
                 correction_method = "bonferroni",
                 domain_scope = "known",
                 evcodes = TRUE,
                 significant = FALSE,
                 sources = c("GO:BP", "REAC", "CORUM", "WP"))$result
  go_res$set <- grp
  enrich_results <- bind_rows(enrich_results, go_res)
}

#### STEP 7: SCORING ENRICHMENT RESULTS ####
enrich_results$FE <- (enrich_results$intersection_size / enrich_results$query_size) /
                     (enrich_results$term_size / enrich_results$effective_domain_size)
enrich_results$stat <- enrich_results$FE * -log10(enrich_results$p_value)
enrich_results <- filter(enrich_results, term_size > 25 & term_size < 500)
enrich_results$genes <- strsplit(as.character(enrich_results$intersection), ",")

#### STEP 8: REDUNDANCY REDUCTION WITH JACCARD FILTER ####
get_top_enrichments <- function(df, set_label, jaccard_threshold = 0.1) {
  df_set <- filter(df, set == set_label)
  df_set <- df_set[order(df_set$p_value), ]
  top_list <- list(df_set[1, ])
  for (i in 2:nrow(df_set)) {
    cur_genes <- df_set$genes[[i]]
    overlaps <- sapply(top_list, function(existing) {
      length(intersect(cur_genes, existing$genes[[1]])) / length(union(cur_genes, existing$genes[[1]]))
    })
    if (all(overlaps < jaccard_threshold)) {
      top_list <- c(top_list, list(df_set[i, ]))
    }
  }
  return(bind_rows(top_list))
}

final_enrichments <- bind_rows(lapply(unique(enrich_results$set), function(s) {
  get_top_enrichments(enrich_results, s)
}))

#### STEP 9: EXPORT ENRICHMENT RESULTS ####
write_csv(final_enrichments, "output/flavi_cluster_enrichment_DIS_filtered.csv")
