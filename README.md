# YFV-human-aedes README

Yellow Fever Virus Interactomes Reveal Common and Divergent Strategies of Replication and Evolution for Mosquito-borne Flaviviruses
  - <https://doi.org/10.1101/2025.06.14.659623>

## Code for data analyses
### MiST
  - 23 JUN 2025 - (Uploaded) humanmist.R and aedesmist.R

    MiST analyses for the YFV-human and YFV-mosquito interactomes. This data uses standard MiST setup files to perform the analysis.
    The MiST package was developed by the Krogan Lab at UCSF: <https://doi.org/10.1002/0471250953.bi0819s49.>

### Enrichment and Clustering
  - 23 JUN 2025 - (Uploaded) enrichment.R and cluster.R

    Enrichment and clustering analyses for Gene Ontology and differential interaction scoring (DIS) for comparative analysis of flavivirus MiST data.

### Foldseek, FATCAT, and NEEDLE analysis
  - 23 JUN 2025 - (Uploaded) Scripts used to generate Foldseek, FATCAT, and NEEDLE data.

    Various scripts used to perform structure and sequence alignments and write to data files.

  - 23 JUN 2025 - (Uploaded) YFV_DENV_Homologs.R

    R script used to evaluate the overlaps between YFV-APMS and DENV-APMS results to identify homolog pairs and interologs.

### Flavivirus conservation analysis
  - 23 JUN 2025 - (Uploaded) Conservation_Analysis_Bootstrapping.R

    Performs correlative analysis of flavivirus proteins and human PPIs
  
## Data for human-Aedes structural and sequence alignments
  - 12 MAR 2025 - (Uploaded) YFV-Foldseekout-Forward.tar.gz and YFV-Foldseek-out-Reverse.tar.gz

    Full .tar.gz Foldseek outputs are provided for the forward (Human to Aedes aegypti) and reverse (Aedes aegypti to Human) searches.
    The search results were limited to e-value < 1e-10.
    Blank text files were generated because a Foldseek search occurred but no results met the e-value threshold.
  
  - ~25 MAR 2025 - (Uploaded) YFV-NEEDLE-out.tar.gz~

    ~Full .tar.gz NEEDLE outputs for the Human and Aedes aegypti proteins.~
  
  - 23 JUN 2025 - (Uploaded) YFV-NEEDLE-out-1.tar.gz and YFV-NEEDLE-out-2.tar.gz

    Replaced and updated previous NEEDLE output archive from 25 MAR 2025. Archive size was too large and was split into two files.

## Data for flavivirus conservation analysis
  - 23 JUN 2025 - (Uploaded) Capsid MiST.csv, NS2B3 MiST.csv, NS3 MiST.csv, NS4A MiST.csv, and NS5 MiST.csv, and Flavivirus_MiST_Correlation.csv

    Full datasets for MiST scores of different flaviviruses for the respective viral proteins, along with correlation data.


Copyright (C) 2025 Shah Lab, UC Davis

E-mail: <prsshah@ucdavis.edu>

Copyright License: GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

