### Human and Aedes YFV and DENV AP-MS Homolog Search
### Chase Skawinski, 05 DEC 2024
### Last Edited, 28 MAY 2025
### Description: Accumulate and define structural homologs using Foldseek, FATCAT, and NEEDLE

### Copyright (C) 2025, Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu


library(dplyr)
library(stringr)
library(ggplot2)
library(ggridges)
library(ggforce)
library(tidyr)

setwd("YOUR_PATH")
getwd()


# Determine Human-Aedes Homologs from Foldseek Search for FATCAT Alignments ------------------------------
  ### Read .csv files
    ## Foldseek data
      # Human to Aedes homologs from Foldseek
        data_foldseek_htoa = read.csv("data/foldseek-Human_to_Aedes-Full-20240823.csv", stringsAsFactors = FALSE)
      # Aedes to Human homologs from Foldseek
        data_foldseek_atoh = read.csv("data/foldseek-Aedes_to_Human-Full-20241007.csv", stringsAsFactors = FALSE)
  
  ### Filter Foldseek hits
    ## Identify overlapping homologs between both searches
      # Clean up the Foldseek information to only contain Uniprot IDs
        data_foldseek_htoa$query = gsub("AF-|-model_v4|\\..*", "", data_foldseek_htoa$query)
        data_foldseek_htoa$target = gsub("AF-|-model_v4|\\..*", "", data_foldseek_htoa$target)
      # Clean up the Foldseek information to only contain Uniprot IDs
        data_foldseek_atoh$query = gsub("AF-|-model_v4|\\..*", "", data_foldseek_atoh$query)
        data_foldseek_atoh$target = gsub("AF-|-model_v4|\\..*", "", data_foldseek_atoh$target)
      # Compare human to Aedes and pull the common homologs
        foldseek_matched = data_foldseek_htoa %>%
          semi_join(data_foldseek_atoh, by = c("query" = "target", "target" = "query")) %>%
          distinct()
      # Add a column of the Uniprot id's
        foldseek_matched$Human.Uniprot = sub("^(.*?)-F\\d.*$", "\\1", foldseek_matched$query)
        foldseek_matched$Aedes.Uniprot = sub("^(.*?)-F\\d.*$", "\\1", foldseek_matched$target)
      # Add a column of the Pair
        foldseek_matched$Pair = paste0(foldseek_matched$Human.Uniprot," ",foldseek_matched$Aedes.Uniprot)
      # Save the data
        write.csv(foldseek_matched, "data/foldseek_matched.csv", row.names = FALSE, quote = FALSE)
    ## Filter the data
      # Filter the human to Aedes file for the 3 lowest e-value matches
        foldseek_matched_filt = foldseek_matched %>%
          group_by(query) %>%
          arrange(evalue) %>%
          slice_head(n = 3)
      # Save the data
        write.csv(foldseek_matched_filt, "data/foldseek_matched_filt.csv", row.names = FALSE, quote = FALSE)


# Combine Foldseek and inParanoid8 FATCAT data -----------------------------
  ## FATCAT homolog searches
    # FATCAT run of all inParanoid8 homologs
      data_fatcat_inParanoid8 = read.csv("data/FATCAT-Human_Aedes-inParanoid8-20241003.csv", stringsAsFactors = FALSE)
      data_fatcat_inParanoid8$Group = NULL
    # FATCAT run of Human-Aedes homologs predicted by Foldseek
      data_fatcat_foldseek = read.csv("data/FATCAT-alignment-foldseek_matched-compfirst.csv", stringsAsFactors = FALSE)
    # Combine the inParanoid-FATCAT and Foldseek-FATCAT data frames (makes things easy to handle)
      data_fatcat = rbind(data_fatcat_inParanoid8, data_fatcat_foldseek)
    # Remove missing data/extraneous rows (i.e., remove inParanoid homologs that don't have structures)
      data_fatcat = data_fatcat[!is.na(data_fatcat$ini.len), ]
    # Remove unnecessary source data frames
      rm(data_fatcat_inParanoid8, data_fatcat_foldseek)
  
  ## Trim the FATCAT protein IDs from the long AF file name
    data_fatcat$Human.Protein = gsub("AF-|-model_v4|\\..*", "", data_fatcat$Human.Protein)
    data_fatcat$Aedes.Protein = gsub("AF-|-model_v4|\\..*", "", data_fatcat$Aedes.Protein)
    # Create a new column with only the Human Uniprot entries
      data_fatcat$Human.Uniprot = sub("-F\\d.*$", "\\1", data_fatcat$Human.Protein)
    # Create a new column with only the Aedes Uniprot entries
      data_fatcat$Aedes.Uniprot = sub("-F\\d.*$", "\\1", data_fatcat$Aedes.Protein)
  ## Make sure there are no duplicated rows from overlapping FATCAT runs
    data_fatcat = unique(data_fatcat)
    # Add a column for Pairs
    data_fatcat$Pair = paste0(data_fatcat$Human.Uniprot," ",data_fatcat$Aedes.Uniprot)
  
  ## Add the source as inParanoid8, Foldseek, or Both
    # Filtered Foldseek source data
      data_foldseek = read.csv("data/foldseek_matched_filt.csv", stringsAsFactors = FALSE)
        data_foldseek$Pair = paste0(data_foldseek$Human.Uniprot," ",data_foldseek$Aedes.Uniprot)
    # inParanoid8 source data
      data_inparanoid8 = read.csv("data/inParanoid8-Human_Aedes-Simple.csv", stringsAsFactors = FALSE)
        data_inparanoid8$Pair = paste0(data_inparanoid8$Protein1," ",data_inparanoid8$Protein2)
    # Create a set of data to test for Source
      both_pairs = intersect(data_foldseek$Pair, data_inparanoid8$Pair)
      foldseek_only = setdiff(data_foldseek$Pair, both_pairs)
      inparanoid_only = setdiff(data_inparanoid8$Pair, both_pairs)
    # Assign Source as one of three options
      data_fatcat = data_fatcat %>%
        mutate(Source = case_when(
            Pair %in% both_pairs ~ "Both",
            Pair %in% foldseek_only ~ "Foldseek",
            Pair %in% inparanoid_only ~ "inParanoid8",
            TRUE ~ NA_character_))
      rm(both_pairs, foldseek_only, inparanoid_only)

  ## Write to a .csv file
    write.csv(data_fatcat, "data/FATCAT-Human_Aedes-inParanoid_Foldseek-Combined.csv", row.names = FALSE, quote = FALSE)


# Group Homolog Results from FATCAT Based on Source (inParanoid8 or Foldseek) ----------------------------------------
  ### Read .csv files
    ## YFV AP-MS Hits Files
      # YFV-Human hits file
        hits_human_yfv = read.csv("data/YFV-Human-Hits.csv", stringsAsFactors = FALSE)
      # YFV-Aedes hits file
        hits_aedes_yfv = read.csv("data/YFV-Aedes-Hits.csv", stringsAsFactors = FALSE)
    ## DENV AP-MS Hits Files
      # DENV-Human hits file
        hits_human_denv = read.csv("data/DENV-Human-Hits.csv", stringsAsFactors = FALSE)
      # DENV-Human hits file
        hits_aedes_denv = read.csv("data/DENV-Aedes-Hits.csv", stringsAsFactors = FALSE)    
    ## Matched Foldseek and inParanoid8 data
      data_fatcat = read.csv("data/FATCAT-Human_Aedes-inParanoid_Foldseek-Combined.csv", stringsAsFactors = FALSE)
    ## inParanoid8 data
      data_inparanoid8 = read.csv("data/inParanoid8-Human_Aedes-Homologs.csv", stringsAsFactors = FALSE)
    ## Matched Foldseek data
      data_foldseek = read.csv("data/foldseek_matched_filt.csv", stringsAsFactors = FALSE)
    
  ### YFV-Human Homologs
    ## Human inParanoid8 Homologs
      # Write inParanoid8 homologs into the hits file
        homologs_human_yfv = hits_human_yfv %>%
          left_join(data_inparanoid8, by = c("Uniprot" = "Protein.ID"))
        homologs_human_yfv = homologs_human_yfv %>%
          filter(!is.na(Group)) %>%
          group_by(Uniprot, Group, Organism) %>%
          summarise(inParanoid8_Homologs = paste(setdiff(data_inparanoid8$Protein.ID[
            data_inparanoid8$Group == Group & data_inparanoid8$Organism != Organism], Uniprot), collapse = ";"), .groups = 'drop')
      # Combine the homologs back into the hits data frame
        hits_human_yfv = hits_human_yfv %>%
          left_join(homologs_human_yfv, by = "Uniprot")
        rm(homologs_human_yfv)
      # Remove the organism column that got added in
        hits_human_yfv$Organism = NULL
      # Change the column name for the group to distinguish as inParanoid8
        hits_human_yfv = hits_human_yfv %>% rename(inParanoid8_Group = Group)
    ## Human Foldseek Homologs (Consensus across Foldseek of Human-->Aedes and Aedes-->Human)
      # Write Foldseek homologs into the hits file
        homologs_human_yfv = hits_human_yfv %>%
          left_join(data_foldseek, by = c("Uniprot" = "Human.Uniprot")) %>%
          group_by(Uniprot) %>%
          summarise(Foldseek_Homologs_Match = paste(na.omit(unique(Aedes.Uniprot)), collapse = ";"), .groups = 'drop')
    ## Combine the homologs back into the hits data frame
      hits_human_yfv = hits_human_yfv %>%
        left_join(homologs_human_yfv, by = "Uniprot")
      rm(homologs_human_yfv)
    ## Write out overlapping inParanoid8 and Foldseek homologs
      hits_human_yfv = hits_human_yfv %>%
        rowwise() %>%
        mutate(
          Both_Homologs = list(intersect(str_split(inParanoid8_Homologs, ";")[[1]], str_split(Foldseek_Homologs_Match, ";")[[1]])),
          Both_Homologs = paste(Both_Homologs, collapse = ";")) %>%
        ungroup()
    ## Write to a .csv file
      write.csv(hits_human_yfv, "data/YFV-Human-Hits_Matched.csv", row.names = FALSE, quote = FALSE)
    
  #### YFV-Aedes Homologs 
    ### Aedes inParanoid8 Homologs
      ## Write inParanoid8 homologs into the hits file for Prey
        homologs_aedes_yfv = hits_aedes_yfv %>%
          left_join(data_inparanoid8, by = c("Prey" = "Protein.ID"))
        homologs_aedes_yfv = homologs_aedes_yfv %>%
          filter(!is.na(Group)) %>%
          group_by(Prey, Group, Organism) %>%
          summarise(inParanoid8_Homologs_Prey = paste(setdiff(data_inparanoid8$Protein.ID[
            data_inparanoid8$Group == Group & data_inparanoid8$Organism != Organism], Prey), collapse = ";"), .groups = 'drop')
        # Combine the homologs back into the hits data frame
          hits_aedes_yfv = hits_aedes_yfv %>%
            left_join(homologs_aedes_yfv, by = "Prey")
          rm(homologs_aedes_yfv)
        # Remove the organism column that got added in
          hits_aedes_yfv$Organism = NULL
        # Change the column name for the group to distinguish as inParanoid8
          hits_aedes_yfv = hits_aedes_yfv %>% rename(inParanoid8_Group_Prey = Group)
      ## Write inParanoid8 homologs into the hits file for PreyGene
        homologs_aedes_yfv = hits_aedes_yfv %>%
          left_join(data_inparanoid8, by = c("PreyGene" = "Protein.ID"))
        homologs_aedes_yfv = homologs_aedes_yfv %>%
          filter(!is.na(Group)) %>%
          group_by(PreyGene, Group, Organism) %>%
          summarise(inParanoid8_Homologs_PreyGene = paste(setdiff(data_inparanoid8$Protein.ID[
            data_inparanoid8$Group == Group & data_inparanoid8$Organism != Organism], PreyGene), collapse = ";"), .groups = 'drop')
        # Combine the homologs back into the hits data frame
          hits_aedes_yfv = hits_aedes_yfv %>%
            left_join(homologs_aedes_yfv, by = "PreyGene")
          rm(homologs_aedes_yfv)
        # Remove the organism column that got added in
          hits_aedes_yfv$Organism = NULL
        # Change the column name for the group to distinguish as inParanoid8
          hits_aedes_yfv = hits_aedes_yfv %>% rename(inParanoid8_Group_PreyGene = Group)
      ## Generate a consensus column in order to have all matches if desired
        hits_aedes_yfv = hits_aedes_yfv %>%
          rowwise() %>%
          mutate(inParanoid8_Homologs = paste(unique(na.omit(c(inParanoid8_Homologs_Prey, inParanoid8_Homologs_PreyGene))), collapse = ";"),
                 inParanoid8_Homologs = str_remove_all(inParanoid8_Homologs, "^;|;$")) %>%
          ungroup()
    ### Aedes Foldseek Homologs (Consensus across Foldseek of Human-->Aedes and Aedes-->Human)
    ### Issue is that Prey and Preygene columns are not the same. We'll make separate columns and a consensus column
      ## Find matches in the Prey column
        data_aedes_prey = hits_aedes_yfv %>%
          left_join(data_foldseek, by = c("Prey" = "Aedes.Uniprot")) %>%
          group_by(Prey) %>%
          summarise(Foldseek_Homologs_Prey = paste(na.omit(unique(Human.Uniprot)), collapse = ";"), .groups = 'drop')
        hits_aedes_yfv = hits_aedes_yfv %>%
          left_join(data_aedes_prey, by = "Prey")
        rm(data_aedes_prey)
      ## Find matches in the PreyGene column
        data_aedes_preygene = hits_aedes_yfv %>%
          left_join(data_foldseek, by = c("PreyGene" = "Aedes.Uniprot")) %>%
          group_by(PreyGene) %>%
          summarise(Foldseek_Homologs_PreyGene = paste(na.omit(unique(Human.Uniprot)), collapse = ";"), .groups = 'drop')
        hits_aedes_yfv = hits_aedes_yfv %>%
          left_join(data_aedes_preygene, by = "PreyGene")
        rm(data_aedes_preygene)
      ## Generate a consensus column in order to have all matches if desired
        hits_aedes_yfv = hits_aedes_yfv %>%
          rowwise() %>%
          mutate(Foldseek_Homologs = paste(unique(na.omit(c(Foldseek_Homologs_Prey, Foldseek_Homologs_PreyGene))), collapse = ";"),
                 Foldseek_Homologs = str_remove_all(Foldseek_Homologs, "^;|;$")) %>%
          ungroup()
      ## Write out overlapping inParanoid8 and Foldseek homologs
        hits_aedes_yfv = hits_aedes_yfv %>%
          rowwise() %>%
          mutate(
            Both_Homologs = list(intersect(str_split(inParanoid8_Homologs, ";")[[1]], str_split(Foldseek_Homologs, ";")[[1]])),
            Both_Homologs = paste(Both_Homologs, collapse = ";")) %>%
          ungroup()
    ### Write to a .csv file
      write.csv(hits_aedes_yfv, "data/YFV-Aedes-Hits_Matched.csv", row.names = FALSE, quote = FALSE)
    
    ### DENV-Human Homologs
      ## Human inParanoid8 Homologs
        # Write inParanoid8 homologs into the hits file
          homologs_human_denv = hits_human_denv %>%
            left_join(data_inparanoid8, by = c("Uniprot" = "Protein.ID"))
          homologs_human_denv = homologs_human_denv %>%
            filter(!is.na(Group)) %>%
            group_by(Uniprot, Group, Organism) %>%
            summarise(inParanoid8_Homologs = paste(setdiff(data_inparanoid8$Protein.ID[
              data_inparanoid8$Group == Group & data_inparanoid8$Organism != Organism], Uniprot), collapse = ";"), .groups = 'drop')
        # Combine the homologs back into the hits data frame
          hits_human_denv = hits_human_denv %>%
            left_join(homologs_human_denv, by = "Uniprot")
          rm(homologs_human_denv)
        # Remove the organism column that got added in
          hits_human_denv$Organism = NULL
        # Change the column name for the group to distinguish as inParanoid8
          hits_human_denv = hits_human_denv %>% rename(inParanoid8_Group = Group)
      ## Human Foldseek Homologs (Consensus across Foldseek of Human-->Aedes and Aedes-->Human)
        # Write Foldseek homologs into the hits file
          homologs_human_denv = hits_human_denv %>%
            left_join(data_foldseek, by = c("Uniprot" = "Human.Uniprot")) %>%
            group_by(Uniprot) %>%
            summarise(Foldseek_Homologs_Match = paste(na.omit(unique(Aedes.Uniprot)), collapse = ";"), .groups = 'drop')
      ## Combine the homologs back into the hits data frame
        hits_human_denv = hits_human_denv %>%
          left_join(homologs_human_denv, by = "Uniprot")
        rm(homologs_human_denv)
      ## Write out overlapping inParanoid8 and Foldseek homologs
        hits_human_denv = hits_human_denv %>%
          rowwise() %>%
          mutate(
            Both_Homologs = list(intersect(str_split(inParanoid8_Homologs, ";")[[1]], str_split(Foldseek_Homologs_Match, ";")[[1]])),
            Both_Homologs = paste(Both_Homologs, collapse = ";")) %>%
          ungroup()
      ## Write to a .csv file
        write.csv(hits_human_denv, "data/DENV-Human-Hits_Matched.csv", row.names = FALSE, quote = FALSE)
  
    ### DENV-Aedes Homologs 
      ## Aedes inParanoid8 Homologs
        # Write inParanoid8 homologs into the hits file
          homologs_aedes_denv = hits_aedes_denv %>%
            left_join(data_inparanoid8, by = c("Uniprot" = "Protein.ID"))
          homologs_aedes_denv = homologs_aedes_denv %>%
            filter(!is.na(Group)) %>%
            group_by(Uniprot, Group, Organism) %>%
            summarise(inParanoid8_Homologs = paste(setdiff(data_inparanoid8$Protein.ID[
              data_inparanoid8$Group == Group & data_inparanoid8$Organism != Organism], Uniprot), collapse = ";"), .groups = 'drop')
        # Combine the homologs back into the hits data frame
          hits_aedes_denv = hits_aedes_denv %>%
            left_join(homologs_aedes_denv, by = "Uniprot")
          rm(homologs_aedes_denv)
        # Remove the organism column that got added in
          hits_aedes_denv$Organism = NULL
        # Change the column name for the group to distinguish as inParanoid8
          hits_aedes_denv = hits_aedes_denv %>% rename(inParanoid8_Group = Group)
      ## Aedes Foldseek Homologs (Consensus across Foldseek of Human-->Aedes and Aedes-->Human)
        homologs_aedes_denv = hits_aedes_denv %>%
          left_join(data_foldseek, by = c("Uniprot" = "Aedes.Uniprot")) %>%
          group_by(Uniprot) %>%
          summarise(Foldseek_Homologs_Match = paste(na.omit(unique(Human.Uniprot)), collapse = ";"), .groups = 'drop')
        hits_aedes_denv = hits_aedes_denv %>%
          left_join(homologs_aedes_denv, by = "Uniprot")
        rm(homologs_aedes_denv)
      ## Write out overlapping inParanoid8 and Foldseek homologs
        hits_aedes_denv = hits_aedes_denv %>%
          rowwise() %>%
          mutate(
            Both_Homologs = list(intersect(str_split(inParanoid8_Homologs, ";")[[1]], str_split(Foldseek_Homologs_Match, ";")[[1]])),
            Both_Homologs = paste(Both_Homologs, collapse = ";")) %>%
          ungroup()
      ## Write to a .csv file
        write.csv(hits_aedes_denv, "data/DENV-Aedes-Hits_Matched.csv", row.names = FALSE, quote = FALSE)


# Categorize Interologs Based on FATCAT and AP-MS Results ----------------------------------------------------
  ### Read .csv files
    ## YFV AP-MS Matched Hits Files
      # YFV-Human matched hits file
        hits_human_yfv = read.csv("data/YFV-Human-Hits_Matched.csv", stringsAsFactors = FALSE)
      # YFV-Aedes matched hits file
        hits_aedes_yfv = read.csv("data/YFV-Aedes-Hits_Matched.csv", stringsAsFactors = FALSE)
    ## DENV AP-MS Matched Hits Files
      # DENV-Human matched hits file
        hits_human_denv = read.csv("data/DENV-Human-Hits_Matched.csv", stringsAsFactors = FALSE)
      # DENV-Human matched hits file
        hits_aedes_denv = read.csv("data/DENV-Aedes-Hits_Matched.csv", stringsAsFactors = FALSE)
    ## FATCAT file from inParanoid8 and Foldseek homologs
      data_fatcat = read.csv("data/FATCAT-Human_Aedes-inParanoid_Foldseek-Combined.csv", stringsAsFactors = FALSE)
  
  ### Identify overlapping interactors with the YFV AP-MS datasets
    ## Match Human PPIs with FATCAT
      data_human = data_fatcat %>%
        left_join(hits_human_yfv, by = c("Human.Uniprot" = "Uniprot")) %>%
        group_by(Human.Uniprot) %>%
        summarise(YFV_Human_Hits = paste(na.omit(unique(Bait)), collapse = ";"), .groups = 'drop')
      # Combine the hits back into the FATCAT data frame
        data_fatcat = data_fatcat %>%
          left_join(data_human, by = "Human.Uniprot")
        rm(data_human)
    ## Match Aedes PPIs with FATCAT (Prey and PreyGene columns don't match. We will combine them.)
      # Find matches in the Prey column
        data_aedes_prey = data_fatcat %>%
          left_join(hits_aedes_yfv, by = c("Aedes.Uniprot" = "Prey")) %>%
          group_by(Aedes.Uniprot) %>%
          summarise(YFV_Aedes_Prey_Hits = paste(na.omit(unique(Bait)), collapse = ";"), .groups = 'drop')
        data_fatcat = data_fatcat %>%
          left_join(data_aedes_prey, by = "Aedes.Uniprot")
        rm(data_aedes_prey)
      # Find matches in the PreyGene column
        data_aedes_preygene = data_fatcat %>%
          left_join(hits_aedes_yfv, by = c("Aedes.Uniprot" = "PreyGene")) %>%
          group_by(Aedes.Uniprot) %>%
          summarise(YFV_Aedes_PreyGene_Hits = paste(na.omit(unique(Bait)), collapse = ";"), .groups = 'drop')
        data_fatcat = data_fatcat %>%
          left_join(data_aedes_preygene, by = "Aedes.Uniprot")
        rm(data_aedes_preygene)
      # Generate a consensus column in order to have all matches if desired
        data_fatcat = data_fatcat %>%
          rowwise() %>%
          mutate(YFV_Aedes_Hits = paste(unique(na.omit(c(YFV_Aedes_Prey_Hits, YFV_Aedes_PreyGene_Hits))), collapse = ";"),
                 YFV_Aedes_Hits = str_remove_all(YFV_Aedes_Hits, "^;|;$")) %>%
          ungroup()
      # Remove the Prey and PreyGene columns because we don't need them
        data_fatcat$YFV_Aedes_Prey_Hits = NULL
        data_fatcat$YFV_Aedes_PreyGene_Hits = NULL
  
  ### Identify overlapping interactors with the DENV AP-MS datasets
    ## Match Human PPIs with FATCAT
      data_human = data_fatcat %>%
        left_join(hits_human_denv, by = c("Human.Uniprot" = "Uniprot")) %>%
        group_by(Human.Uniprot) %>%
        summarise(DENV_Human_Hits = paste(na.omit(unique(Bait)), collapse = ";"), .groups = 'drop')
      # Combine the hits back into the FATCAT data frame
        data_fatcat = data_fatcat %>%
          left_join(data_human, by = "Human.Uniprot")
        rm(data_human)
    ## Match Aedes PPIs with FATCAT
      data_aedes = data_fatcat %>%
        left_join(hits_aedes_denv, by = c("Aedes.Uniprot" = "Uniprot")) %>%
        group_by(Aedes.Uniprot) %>%
        summarise(DENV_Aedes_Hits = paste(na.omit(unique(Bait)), collapse = ";"), .groups = 'drop')
      # Combine the hits back into the FATCAT data frame
        data_fatcat = data_fatcat %>%
          left_join(data_aedes, by = "Aedes.Uniprot")
        rm(data_aedes)
  
  
  ###
    ## Eliminate duplicate protein matches (mostly affects proteins with AF structures for multiple sections)
      # Eliminate duplicate pairs by smallest p-value
        data_fatcat_trim = data_fatcat %>%
          group_by(Pair) %>%
          slice_min(order_by = P.value, with_ties = TRUE) %>%
          slice_min(order_by = opt.rmsd, with_ties = FALSE) %>%
          ungroup()
  ###
  
  
  ### Categorize homolog pairs as homologs, HC-PPIs, or interologs
    ## YFV AP-MS data (using combined Prey and PreyGene column)
      data_fatcat_trim = data_fatcat_trim %>%
        mutate(
          # Split the strings into lists of vectors
            vec1 = str_split(YFV_Human_Hits, ";"),
            vec2 = str_split(YFV_Aedes_Hits, ";"),
          # Remove empty strings from the vectors
            vec1 = lapply(vec1, \(x) x[x != ""]),
            vec2 = lapply(vec2, \(x) x[x != ""]),
          # Determine the classification (homolog, HC-PPI, or interolog)
            YFV_Classification = case_when(
              sapply(vec1, length) == 0 & sapply(vec2, length) == 0 ~ "Homolog",
              sapply(vec1, length) > 0 & sapply(vec2, length) > 0 & mapply(function(v1, v2) any(v1 %in% v2), vec1, vec2) ~ "Interolog",
              sapply(vec1, length) > 0 | sapply(vec2, length) > 0 ~ "HC-PPI",
              TRUE ~ NA_character_
            )
          ) %>%
        # Clean up temporary columns
          select(-vec1, -vec2)
    ## YFV AP-MS data (using combined Prey and PreyGene column)
      data_fatcat_trim = data_fatcat_trim %>%
        mutate(
          # Split the strings into lists of vectors
            vec1 = str_split(DENV_Human_Hits, ";"),
            vec2 = str_split(DENV_Aedes_Hits, ";"),
          # Remove empty strings from the vectors
            vec1 = lapply(vec1, \(x) x[x != ""]),
            vec2 = lapply(vec2, \(x) x[x != ""]),
          # Determine the classification (homolog, HC-PPI, or interolog)
            DENV_Classification = case_when(
              sapply(vec1, length) == 0 & sapply(vec2, length) == 0 ~ "Homolog",
              sapply(vec1, length) > 0 & sapply(vec2, length) > 0 & mapply(function(v1, v2) any(v1 %in% v2), vec1, vec2) ~ "Interolog",
              sapply(vec1, length) > 0 | sapply(vec2, length) > 0 ~ "HC-PPI",
              TRUE ~ NA_character_
            )
          ) %>%
        # Clean up temporary columns
          select(-vec1, -vec2)
    
    ## Use median Z-score to standardize RMSD data
      # Add column of Z-scores for ini.rmsd and opt.rmsd
        data_fatcat_trim$z.ini = ((data_fatcat_trim$ini.rmsd - median(data_fatcat_trim$ini.rmsd)) / sd(data_fatcat_trim$ini.rmsd))
        data_fatcat_trim$z.opt = ((data_fatcat_trim$opt.rmsd - median(data_fatcat_trim$opt.rmsd)) / sd(data_fatcat_trim$opt.rmsd))
    
    ## Write the data to a .csv file
      write.csv(data_fatcat_trim, "data/FATCAT-Human_Aedes-Interactors.csv", row.names = FALSE, quote = FALSE)


# Structural Homolog Statistical Significance Tests ----------------------------------------------------------
  ## Read .csv file
    data_wilcoxon = read.csv("data/FATCAT-Human_Aedes-Interactors.csv", stringsAsFactors = FALSE)
    
  ## Perform Wilcoxon Rank-Sum Tests for Initial RMSDs
    # YFV Interologs vs. Other Homologs  
      wilcox_YFV_Interologs_vs_Homologs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification != "Interolog"]
        )
    # DENV Interologs vs. Other Homologs
      wilcox_DENV_Interologs_vs_Homologs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # All Interologs vs. Other Homologs
      wilcox_All_Interologs_vs_Homologs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
        )
    # YFV HC-PPIs vs. Other Homologs
      wilcox_YFV_PPIs_vs_Homologs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$YFV_Classification == "HC-PPI"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification != "Interolog" | data_wilcoxon$YFV_Classification != "HC-PPI"]
      )
    # DENV HC-PPIs vs. Other Homologs
      wilcox_DENV_PPIs_vs_Homologs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$DENV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$DENV_Classification != "Interolog" | data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All HC-PPIs vs. Other Homologs
      wilcox_All_PPIs_vs_Homologs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog" | 
                                 data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog" & 
                                 data_wilcoxon$YFV_Classification != "HC-PPI" & data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All Interologs vs. Other HC-PPIs
      wilcox_All_Ints_vs_PPIs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI" &
                                 data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV Interologs vs. Other YFV HC-PPIs
      wilcox_YFV_Ints_vs_YFV_PPIs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$YFV_Classification == "HC-PPI" & data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other DENV HC-PPIs
      wilcox_DENV_Ints_vs_DENV_PPIs_ini = wilcox.test(
        data_wilcoxon$ini.rmsd[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$ini.rmsd[data_wilcoxon$DENV_Classification == "HC-PPI" & data_wilcoxon$DENV_Classification != "Interolog"]
      )

  ## Perform Wilcoxon Rank-Sum Tests for Optimized RMSDs
    # YFV Interologs vs. Other Homologs  
      wilcox_YFV_Interologs_vs_Homologs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other Homologs
      wilcox_DENV_Interologs_vs_Homologs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # All Interologs vs. Other Homologs
      wilcox_All_Interologs_vs_Homologs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV HC-PPIs vs. Other Homologs
      wilcox_YFV_PPIs_vs_Homologs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$YFV_Classification == "HC-PPI"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification != "Interolog" | data_wilcoxon$YFV_Classification != "HC-PPI"]
      )
    # DENV HC-PPIs vs. Other Homologs
      wilcox_DENV_PPIs_vs_Homologs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$DENV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$DENV_Classification != "Interolog" | data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All HC-PPIs vs. Other Homologs
      wilcox_All_PPIs_vs_Homologs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog" | 
                                 data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog" & 
                                 data_wilcoxon$YFV_Classification != "HC-PPI" & data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All Interologs vs. Other HC-PPIs
      wilcox_All_Ints_vs_PPIs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI" &
                                 data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV Interologs vs. Other YFV HC-PPIs
      wilcox_YFV_Ints_vs_YFV_PPIs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$YFV_Classification == "HC-PPI" & data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other DENV HC-PPIs
      wilcox_DENV_Ints_vs_DENV_PPIs_opt = wilcox.test(
        data_wilcoxon$opt.rmsd[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$opt.rmsd[data_wilcoxon$DENV_Classification == "HC-PPI" & data_wilcoxon$DENV_Classification != "Interolog"]
      )

  ## Combine all results into a single data frame and perform p-value adjustment with Bonferroni correction
    # Create a list of the tests' names
      test_names = ls(pattern = "^wilcox_")
    # Combine the results
      pvalues = do.call(rbind, lapply(test_names, function(name) {
        test_result = get(name)
        data.frame(Test = name, P_Value = test_result$p.value, stringsAsFactors = FALSE)
      }))
      rm(list = ls(pattern = "^wilcox_"))
    # Clean up test names and separate initial from optimized RMSDs
      pvalues$Base_Test = sub("_ini|_opt", "", pvalues$Test)
      pvalues$Base_Test = sub("wilcox_", "", pvalues$Base_Test)
      pvalues = merge(
        pvalues[grep("_ini$", pvalues$Test),],
        pvalues[grep("_opt$", pvalues$Test),],
        by = "Base_Test",
        suffixes = c("_ini", "_opt")
      )
      pvalues$Test_ini = NULL
      pvalues$Test_opt = NULL
    # Adjust p-values with Bonferroni correction
      pvalues$AdjP_Value_ini = p.adjust(pvalues$P_Value_ini, method = "bonferroni")
      pvalues$AdjP_Value_opt = p.adjust(pvalues$P_Value_opt, method = "bonferroni")

  ## Write to a .csv file
    write.csv(pvalues, "data/FATCAT-Human_Aedes-Wilcoxon.csv", row.names = FALSE, quote = FALSE)


# Make RMSD Stats for Homologs, HC-PPIs, and Interologs -------------------------------------------------------------------
  ## Read .csv files
    data = read.csv("data/FATCAT-Human_Aedes-Interactors.csv", stringsAsFactors = FALSE)
    data_sig = read.csv("data/FATCAT-Human_Aedes-Wilcoxon.csv", stringsAsFactors = FALSE)
  
  ### Create comprehensive data frame with all of the information in two columns to plot with ggplot2
    ## Create new column for combined groupings
      data_plot = data
      data_plot$Group = "All Homologs"
    ## Extend the data frame to include all different forms of groups
      data_plot = bind_rows(
        # Keep all Homologs
          data_plot %>% mutate(Group = "All Homologs"),
        # Add YFV HC-PPIs
          data_plot %>% filter(YFV_Classification %in% c("HC-PPI", "Interolog")) %>%
            mutate(Group = "YFV HC-PPIs"),
        # Add YFV Interologs group
          data_plot %>% filter(YFV_Classification == "Interolog") %>%
            mutate(Group = "YFV Interologs"),
        # Add DENV HC-PPIs
          data_plot %>% filter(DENV_Classification %in% c("HC-PPI", "Interolog")) %>%
            mutate(Group = "DENV HC-PPIs"),
        # Add DENV Interologs
          data_plot %>% filter(DENV_Classification == "Interolog") %>%
            mutate(Group = "DENV Interologs"),
        # Combine all HC-PPIs
          data_plot %>% filter(YFV_Classification == "Interolog" | YFV_Classification == "HC-PPI" |
                               DENV_Classification == "Interolog" | DENV_Classification == "HC-PPI") %>%
          mutate(Group = "All HC-PPIs"),
        # Combine all Interologs
          data_plot %>% filter(YFV_Classification == "Interolog" | DENV_Classification == "Interolog") %>%
            mutate(Group = "All Interologs")
        )
    ## Set factors that dictates the order in which our categories will plot
      data_plot$Group = factor(data_plot$Group, levels = c("All Homologs", "All HC-PPIs",
                                                           "YFV HC-PPIs", "DENV HC-PPIs",
                                                           "All Interologs",
                                                           "YFV Interologs", "DENV Interologs")
                               )
    ## Create a data frame that includes number of samples, medians, and means for all classifications
      stats_plot = data_plot %>%
        group_by(Group) %>%
        summarize(
          # Quantities
            n = n(),
          # Medians
            median_ini = median(ini.rmsd, na.rm = TRUE),
            median_opt = median(opt.rmsd, na.rm = TRUE),
          # Means
            mean_ini = mean(ini.rmsd, na.rm = TRUE),
            mean_opt = mean(opt.rmsd, na.rm = TRUE),
        )
      # Write to a .csv file
        write.csv(stats_plot, "data/FATCAT-Human_Aedes-Stats.csv", row.names = FALSE, quote = FALSE)


# Combine structural homolog data with NEEDLE alignments -----------------------------
  # Read .csv files
    data = read.csv("data/FATCAT-Human_Aedes-Interactors.csv", stringsAsFactors = FALSE)
    data_needle = read.csv("data/NEEDLE-Human_Aedes-homologs.csv", stringsAsFactors = FALSE)
    
  # Combine Homolog data with NEEDLE data
    data_comb = data %>%
      left_join(data_needle, by = c("Human.Uniprot" = "Human.Protein", "Aedes.Uniprot" = "Aedes.Protein"))
  
  # Save as a new .csv file
    write.csv(data_comb, "data/FATCAT_NEEDLE-Human_Aedes-Interactors.csv", row.names = FALSE, quote = FALSE)


# NEEDLE Alignment Statistical Significance Tests -------------------------
  ## Read .csv file
    data_wilcoxon = read.csv("data/FATCAT_NEEDLE-Human_Aedes-Interactors.csv", stringsAsFactors = FALSE)
    
    data_wilcoxon = data_wilcoxon %>%
      mutate(
        Identity = as.numeric(sub("%", "", Identity)),
        Similarity = as.numeric(sub("%", "", Similarity)),
        Gaps = as.numeric(sub("%", "", Gaps))
      )
  ## Perform Wilcoxon Rank-Sum Tests for AA Identity
    # YFV Interologs vs. Other Homologs
      wilcox_YFV_Interologs_vs_Homologs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other Homologs
      wilcox_DENV_Interologs_vs_Homologs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Identity[data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # All Interologs vs. Other Homologs
      wilcox_All_Interologs_vs_Homologs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV HC-PPIs vs. Other Homologs
      wilcox_YFV_PPIs_vs_Homologs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$YFV_Classification == "HC-PPI"],
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification != "Interolog" | data_wilcoxon$YFV_Classification != "HC-PPI"]
      )
    # DENV HC-PPIs vs. Other Homologs
      wilcox_DENV_PPIs_vs_Homologs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$DENV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$Identity[data_wilcoxon$DENV_Classification != "Interolog" | data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All HC-PPIs vs. Other Homologs
      wilcox_All_PPIs_vs_Homologs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog" | 
                                 data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog" & 
                                 data_wilcoxon$YFV_Classification != "HC-PPI" & data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All Interologs vs. Other HC-PPIs
      wilcox_All_Ints_vs_PPIs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI" &
                                 data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV Interologs vs. Other YFV HC-PPIs
      wilcox_YFV_Ints_vs_YFV_PPIs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$Identity[data_wilcoxon$YFV_Classification == "HC-PPI" & data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other DENV HC-PPIs
      wilcox_DENV_Ints_vs_DENV_PPIs_ident = wilcox.test(
        data_wilcoxon$Identity[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Identity[data_wilcoxon$DENV_Classification == "HC-PPI" & data_wilcoxon$DENV_Classification != "Interolog"]
      )

  ## Perform Wilcoxon Rank-Sum Tests for AA Similarity
    # YFV Interologs vs. Other Homologs  
      wilcox_YFV_Interologs_vs_Homologs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other Homologs
      wilcox_DENV_Interologs_vs_Homologs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Similarity[data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # All Interologs vs. Other Homologs
      wilcox_All_Interologs_vs_Homologs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV HC-PPIs vs. Other Homologs
      wilcox_YFV_PPIs_vs_Homologs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$YFV_Classification == "HC-PPI"],
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification != "Interolog" | data_wilcoxon$YFV_Classification != "HC-PPI"]
      )
    # DENV HC-PPIs vs. Other Homologs
      wilcox_DENV_PPIs_vs_Homologs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$DENV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$Similarity[data_wilcoxon$DENV_Classification != "Interolog" | data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All HC-PPIs vs. Other Homologs
      wilcox_All_PPIs_vs_Homologs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog" | 
                                 data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog" & 
                                 data_wilcoxon$YFV_Classification != "HC-PPI" & data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All Interologs vs. Other HC-PPIs
      wilcox_All_Ints_vs_PPIs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI" &
                                 data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV Interologs vs. Other YFV HC-PPIs
      wilcox_YFV_Ints_vs_YFV_PPIs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$Similarity[data_wilcoxon$YFV_Classification == "HC-PPI" & data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other DENV HC-PPIs
      wilcox_DENV_Ints_vs_DENV_PPIs_sim = wilcox.test(
        data_wilcoxon$Similarity[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Similarity[data_wilcoxon$DENV_Classification == "HC-PPI" & data_wilcoxon$DENV_Classification != "Interolog"]
      )

  ## Perform Wilcoxon Rank-Sum Tests for AA Gaps
    # YFV Interologs vs. Other Homologs  
      wilcox_YFV_Interologs_vs_Homologs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other Homologs
      wilcox_DENV_Interologs_vs_Homologs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Gaps[data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # All Interologs vs. Other Homologs
      wilcox_All_Interologs_vs_Homologs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV HC-PPIs vs. Other Homologs
      wilcox_YFV_PPIs_vs_Homologs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$YFV_Classification == "HC-PPI"],
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification != "Interolog" | data_wilcoxon$YFV_Classification != "HC-PPI"]
      )
    # DENV HC-PPIs vs. Other Homologs
      wilcox_DENV_PPIs_vs_Homologs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$DENV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$Gaps[data_wilcoxon$DENV_Classification != "Interolog" | data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All HC-PPIs vs. Other Homologs
      wilcox_All_PPIs_vs_Homologs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog" | 
                                   data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI"],
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog" & 
                                   data_wilcoxon$YFV_Classification != "HC-PPI" & data_wilcoxon$DENV_Classification != "HC-PPI"]
      )
    # All Interologs vs. Other HC-PPIs
      wilcox_All_Ints_vs_PPIs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification == "Interolog" | data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification == "HC-PPI" | data_wilcoxon$DENV_Classification == "HC-PPI" &
                                   data_wilcoxon$YFV_Classification != "Interolog" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
    # YFV Interologs vs. Other YFV HC-PPIs
      wilcox_YFV_Ints_vs_YFV_PPIs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification == "Interolog"],
        data_wilcoxon$Gaps[data_wilcoxon$YFV_Classification == "HC-PPI" & data_wilcoxon$YFV_Classification != "Interolog"]
      )
    # DENV Interologs vs. Other DENV HC-PPIs
      wilcox_DENV_Ints_vs_DENV_PPIs_gap = wilcox.test(
        data_wilcoxon$Gaps[data_wilcoxon$DENV_Classification == "Interolog"],
        data_wilcoxon$Gaps[data_wilcoxon$DENV_Classification == "HC-PPI" & data_wilcoxon$DENV_Classification != "Interolog"]
      )
  
  ## Combine all results into a single data frame and perform p-value adjustment with Bonferroni correction
    # Create a list of the tests' names
      test_names = ls(pattern = "^wilcox_")
    # Combine the results
      pvalues = do.call(rbind, lapply(test_names, function(name) {
        test_result = get(name)
        data.frame(Test = name, P_Value = test_result$p.value, stringsAsFactors = FALSE)
      }))
      rm(list = ls(pattern = "^wilcox_"))
    # Clean up test names and separate Identity, Similarity, and Gaps
      pvalues$Base_Test = sub("_ident|_sim|_gap", "", pvalues$Test)
      pvalues$Base_Test = sub("wilcox_", "", pvalues$Base_Test)
      pvalues_ident <- pvalues[grep("_ident$", pvalues$Test),]
      pvalues_sim <- pvalues[grep("_sim$", pvalues$Test),]
      pvalues_gap <- pvalues[grep("_gap$", pvalues$Test),]
      
      pvalues_merged <- merge(pvalues_ident,
                              pvalues_sim,
                              by = "Base_Test",
                              suffixes = c("_ident", "_sim")
                              )
      pvalues_merged <- merge(pvalues_merged,
                              pvalues_gap,
                              by = "Base_Test",
                              suffixes = c("", "_gap")
                              )
      
      pvalues_merged = pvalues_merged %>%
        rename(Test_gap = Test) %>%
        rename(P_Value_gap = P_Value)
      
      pvalues_merged$Test_ident = NULL
      pvalues_merged$Test_sim = NULL
      pvalues_merged$Test_gap = NULL
    # Adjust p-values with Bonferroni correction
      pvalues_merged$AdjP_Value_ident = p.adjust(pvalues_merged$P_Value_ident, method = "bonferroni")
      pvalues_merged$AdjP_Value_sim = p.adjust(pvalues_merged$P_Value_sim, method = "bonferroni")
      pvalues_merged$AdjP_Value_gap = p.adjust(pvalues_merged$P_Value_gap, method = "bonferroni")
  
  ## Write to a .csv file
    write.csv(pvalues_merged, "data/NEEDLE-Human_Aedes-Wilcoxon.csv", row.names = FALSE, quote = FALSE)
  
  
# Make NEEDLE Boxplots for Homologs, HC-PPIs, and Interologs --------------
  ## Read .csv file
    data = read.csv("data/FATCAT_NEEDLE-Human_Aedes-Interactors.csv", stringsAsFactors = FALSE)

  ### Create comprehensive data frame with all of the information in two columns to plot with ggplot2
    ## Create new column for combined groupings
      data_plot = data
      data_plot$Group = "All Homologs"
    ## Extend the data frame to include all different forms of groups
      data_plot = bind_rows(
        # Keep all Homologs
        data_plot %>% mutate(Group = "All Homologs"),
        # Add YFV HC-PPIs
        data_plot %>% filter(YFV_Classification %in% c("HC-PPI", "Interolog")) %>%
          mutate(Group = "YFV HC-PPIs"),
        # Add YFV Interologs group
        data_plot %>% filter(YFV_Classification == "Interolog") %>%
          mutate(Group = "YFV Interologs"),
        # Add DENV HC-PPIs
        data_plot %>% filter(DENV_Classification %in% c("HC-PPI", "Interolog")) %>%
          mutate(Group = "DENV HC-PPIs"),
        # Add DENV Interologs
        data_plot %>% filter(DENV_Classification == "Interolog") %>%
          mutate(Group = "DENV Interologs"),
        # Combine all HC-PPIs
        data_plot %>% filter(YFV_Classification == "Interolog" | YFV_Classification == "HC-PPI" |
                               DENV_Classification == "Interolog" | DENV_Classification == "HC-PPI") %>%
          mutate(Group = "All HC-PPIs"),
        # Combine all Interologs
        data_plot %>% filter(YFV_Classification == "Interolog" | DENV_Classification == "Interolog") %>%
          mutate(Group = "All Interologs")
      )
    ## Set factors that dictates the order in which our categories will plot
      data_plot$Group = factor(data_plot$Group, levels = c("All Homologs", "All HC-PPIs",
                                                           "YFV HC-PPIs", "DENV HC-PPIs",
                                                           "All Interologs",
                                                           "YFV Interologs", "DENV Interologs")
      )
    ## Create a data frame that includes number of samples, medians, and means for all classifications
      # Read the values as numbers, not strings
        data_plot = data_plot %>%
          mutate(
            Identity = as.numeric(sub("%", "", Identity)),
            Similarity = as.numeric(sub("%", "", Similarity)),
            Gaps = as.numeric(sub("%", "", Gaps))
          )
        
      stats_plot = data_plot %>%
        group_by(Group) %>%
        summarize(
          # Quantities
          n = n(),
          # Medians
          median_id = median(Identity, na.rm = TRUE),
          median_sim = median(Similarity, na.rm = TRUE),
          median_gap = median(Gaps, na.rm = TRUE),
          # Means
          mean_id = mean(Identity, na.rm = TRUE),
          mean_sim = mean(Similarity, na.rm = TRUE),
          mean_gap = mean(Gaps, na.rm = TRUE)
        )
    # Write to a .csv file
    write.csv(stats_plot, "data/NEEDLE-Human_Aedes-Stats.csv", row.names = FALSE, quote = FALSE)


# Make RMSD and NEEDLE Boxplots for Pan-Host Clusters -------------------------------------------------------------
  ## Read .csv files
    data = read.csv("data/FATCAT_NEEDLE-Human_Aedes-Interactors.csv", stringsAsFactors = FALSE)
    data_panhost = read.csv("data/data/panhost-clust.csv", stringsAsFactors = FALSE)
  
  ## Create a data frame that has clustered proteins with the respective initial and optimized RMSDs
    # Convert percentages in the NEEDLE data to numbers
      data = data %>%
        mutate(
          Identity = as.numeric(sub("%", "", Identity)),
          Similarity = as.numeric(sub("%", "", Similarity)),
          Gaps = as.numeric(sub("%", "", Gaps))
        )

    # Combine the data into the pan-host clustering data frame
      data_clust = data_panhost %>%
        left_join(data %>% dplyr::select(Human.Uniprot, Aedes.Uniprot, opt.rmsd, Identity),
                  by = c("UNIPROT" = "Human.Uniprot")
                  )

  ## Set factors that dictates the order in which our categories will plot
    data_clust$hclust = factor(data_clust$hclust, levels = c("1", "2", "3", "4",
                                                             "5", "6", "7", "8")
                               )

  ## Create a data frame that includes number of samples, medians, and means for all classifications
    stats_clust = data_clust %>%
      group_by(hclust) %>%
      summarize(
        # Quantities
        n = n(),
        # Medians
        median_opt = median(opt.rmsd, na.rm = TRUE),
        median_ident = median(Identity, na.rm = TRUE),
        # Means
        mean_opt = mean(opt.rmsd, na.rm = TRUE),
        mean_ident = mean(Identity, na.rm = TRUE)
        )

    # Write to a .csv file
      write.csv(stats_clust, "YFV_DENV-PanHost_Clusters-Stats.csv", row.names = FALSE, quote = FALSE)


# Add Panhost cluster 5 to NEEDLE and FATCAT plots ------------------------
  ## Read .csv files
    data = read.csv("data/FATCAT_NEEDLE-Human_Aedes-Interactors.csv", stringsAsFactors = FALSE)
    data_panhost = read.csv("data/data/panhost-clust.csv", stringsAsFactors = FALSE)

  ## Read percentages as numbers
    data = data %>%
      mutate(
        Identity = as.numeric(sub("%", "", Identity)),
        Similarity = as.numeric(sub("%", "", Similarity)),
        Gaps = as.numeric(sub("%", "", Gaps))
      )

  ## Annotate the data for cluster 5
    data_plot = data %>%
      left_join(data_panhost %>% 
                  filter(hclust == 5) %>%
                  dplyr::select(UNIPROT, hclust), by = c("Human.Uniprot" = "UNIPROT")) %>%
      distinct()

  ## Extend the data frame to include all different forms of groups
    data_plot$Group = "All Homologs"
    data_plot = bind_rows(
      # Keep all Homologs
        data_plot %>% mutate(Group = "All Homologs"),
      # Add YFV HC-PPIs
        data_plot %>% filter(YFV_Classification %in% c("HC-PPI", "Interolog")) %>%
          mutate(Group = "YFV HC-PPIs"),
      # Add YFV Interologs group
        data_plot %>% filter(YFV_Classification == "Interolog") %>%
          mutate(Group = "YFV Interologs"),
      # Add DENV HC-PPIs
        data_plot %>% filter(DENV_Classification %in% c("HC-PPI", "Interolog")) %>%
          mutate(Group = "DENV HC-PPIs"),
      # Add DENV Interologs
        data_plot %>% filter(DENV_Classification == "Interolog") %>%
          mutate(Group = "DENV Interologs"),
      # Combine all HC-PPIs
        data_plot %>% filter(YFV_Classification == "Interolog" | YFV_Classification == "HC-PPI" |
                               DENV_Classification == "Interolog" | DENV_Classification == "HC-PPI") %>%
          mutate(Group = "All HC-PPIs"),
      # Combine all Interologs
        data_plot %>% filter(YFV_Classification == "Interolog" | DENV_Classification == "Interolog") %>%
          mutate(Group = "All Interologs"),
      # Use only Pan-host cluster 5
        data_plot %>% filter(hclust == 5) %>%
          mutate(Group = "Pan-Host")
      )
    
  ## Set factors that dictates the order in which our categories will plot
    data_plot$Group = factor(data_plot$Group, levels = c("All Homologs",
                                                         "YFV HC-PPIs", "All HC-PPIs",
                                                         "YFV Interologs", "All Interologs",
                                                         "Pan-Host",
                                                         "DENV HC-PPIs", "DENV Interologs")
                             )

  ### Generate statistical data
    stats_plot = data_plot %>%
      group_by(Group) %>%
      summarize(
        # Quantities
          n = n(),
        ## NEEDLE
          # Sequence Identity
            median_id = median(Identity, na.rm = TRUE),
            mean_id = mean(Identity, na.rm = TRUE),
          # Sequence Similarity
            median_sim = median(Similarity, na.rm = TRUE),
            mean_sim = mean(Similarity, na.rm = TRUE),
          # Sequence Gaps
            median_gap = median(Gaps, na.rm = TRUE),
            mean_gap = mean(Gaps, na.rm = TRUE),
        ## FATCAT
          # Optimized RMSD
            median_opt = median(opt.rmsd, na.rm = TRUE),
            mean_opt = mean(opt.rmsd, na.rm = TRUE),
          # Initial RMSD
            median_ini = median(ini.rmsd, na.rm = TRUE),
            mean_ini = mean(ini.rmsd, na.rm = TRUE)
      )

    # Write to a .csv file
      write.csv(stats_plot, "YFV_DENV-Homolog_PanHost-Stats.csv", row.names = FALSE, quote = FALSE)

  ## Pairwise Wilcoxon tests
    # RMSD, initial and optimized
      adj_p_val_ini.rmsd = pairwise.wilcox.test(data_plot$ini.rmsd, data_plot$Group, p.adjust.method = "bonferroni")
      adj_p_val_opt.rmsd = pairwise.wilcox.test(data_plot$opt.rmsd, data_plot$Group, p.adjust.method = "bonferroni")
    # Amino acid identity, similarity, and gaps
      adj_p_val_identity = pairwise.wilcox.test(data_plot$Identity, data_plot$Group, p.adjust.method = "bonferroni")
      adj_p_val_similarity = pairwise.wilcox.test(data_plot$Similarity, data_plot$Group, p.adjust.method = "bonferroni")
      adj_p_val_gaps = pairwise.wilcox.test(data_plot$Gaps, data_plot$Group, p.adjust.method = "bonferroni")
    # Write adjusted p-values to a set
      test_names = ls(pattern = "^adj_p_val_")
      pval_list = list()
      
      for (name in test_names) {
        test_obj = get(name)
        P_val_adj = test_obj$p.value
        
        P_val_df = as.data.frame(as.table(P_val_adj)) %>%
          filter(!is.na(Freq)) %>%
          rename(Group1 = Var1,
                 Group2 = Var2,
                 !!name := Freq) %>%
          mutate(
            !!paste0(name, "_stars") := case_when(
              !!sym(name) < 0.0001 ~ "****",
              !!sym(name) < 0.001  ~ "***",
              !!sym(name) < 0.01   ~ "**",
              !!sym(name) < 0.05   ~ "*",
              TRUE                 ~ ""
            )
          )
        
        pval_list[[name]] = P_val_df
        }

      P_value_wide <- Reduce(function(x, y) full_join(x, y, by = c("Group1", "Group2")), pval_list)
      
    # Write to a .csv file
      write.csv(P_value_wide, "FATCAT_NEEDLE-Human_Aedes-Wilcoxon.csv", row.names = FALSE, quote = FALSE)
