### Read NEEDLE results and write to a .csv output file
### Chase Skawinski - 23 JUN 2025
### Copyright (C) 2025 Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu

#!/bin/bash

# Set working directory
  echo "Running on host $(hostname)"
  echo "Start time: $(date)"
  echo "Directory: $(pwd)"

# Input path
  needle_dir="needle_out"

# Output .csv file
  output_file="Human_Aedes-homologs-needle.csv"
  # Write the headers for the file
    echo "Human.Protein,Aedes.Protein,Ali.Length,Identity,Similarity,Gaps,Score" > "$output_file"

# Process each .txt file in the input directory
  for file in "$needle_dir"/*.txt; do
    
  # Get human and Aedes protein names
    human_protein=$(grep -m 1 "^# 1: " "$file" | awk '{print $3}')
    aedes_protein=$(grep -m 1 "^# 2: " "$file" | awk '{print $3}')
  # Get alignment length
    ali_length=$(grep -m 1 "^# Length: " "$file" | awk -F': ' '{print $2}' | xargs)
  # Get Identity percentage
    identity=$(grep -m 1 "^# Identity: " "$file" | grep -oP '(?<=\()[^\)]+')
  # Get Similarity percentage
    similarity=$(grep -m 1 "^# Similarity: " "$file" | grep -oP '(?<=\()[^\)]+')
  # Get Gaps percentage
    gaps=$(grep -m 1 "^# Gaps: " "$file" | grep -oP '(?<=\()[^\)]+')
  # Get Score
    score=$(grep -m 1 "^# Score: " "$file" | awk -F': ' '{print $2}' | xargs)
  # Write data to the output file
    echo "$human_protein,$aedes_protein,$ali_length,$identity,$similarity,$gaps,$score" >> "$output_file"
  
  done

echo "Writing of NEEDLE files into ${output_file} is complete"
echo "End time: $(date)"
