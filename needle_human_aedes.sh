### NEEDLE comparison of list of proteins using server
### Reads a .csv file and aligns each protein pair
### Chase Skawinski - 23 JUN 2025
### Copyright (C) 2025 Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu

#!/bin/bash

echo "Running on host $(hostname)"
echo "Start time: $(date)"
echo "Directory: $(pwd)"

# Input .csv file
  input_file="Human_Aedes-Prot_Pairs.csv"

# Output directories
  output_dir="needle_out"
  mkdir -p "$output_dir"
  fasta_dir="fasta"
  mkdir -p "$fasta_dir"

# Loop through .csv file, skip header
  {
  read
  while IFS=, read -r human_protein human_sequence aedes_protein aedes_sequence; do
    
    # Create FASTA files
      human_fasta="${fasta_dir}/${human_protein}.fasta"
      aedes_fasta="${fasta_dir}/${aedes_protein}.fasta"
      
    # Write FASTA sequences to a temp FASTA file
      echo ">${human_protein}" > "$human_fasta"
      echo "$human_sequence" >> "$human_fasta"
      
      echo ">${aedes_protein}" > "$aedes_fasta"
      echo "$aedes_sequence" >> "$aedes_fasta"
      
    # Output alignment file name
      alignment_file="${output_dir}/${human_protein}_vs_${aedes_protein}-needle.txt"
      
    # Run NEEDLE alignment
      needle -asequence "$human_fasta" -bsequence "$aedes_fasta" -gapopen 10.0 -gapextend 0.5 -outfile "$alignment_file"
  done
  } < "$input_file"

echo "Alignments are complete"
echo "End time: $(date)"
