### FATCAT alignments of viral proteins
### Reads a .csv file and aligns each protein against the other
### Change the paths and directories for the comparison needed
### Chase Skawinski - 23 JUN 2025
### Copyright (C) 2025 Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu

#!/bin/bash

# Protein file
  viral_prots="PATH_TO_FILE"
  # Ensure that the file ends with a newline if missing
    if [ -s "${viral_prots}" ] && [ "$(tail -c1 "${viral_prots}")" != "" ]; then
      echo "" >> "${viral_prots}"
    fi
  # Generate arrays for viruses and proteins
    declare -a viruses
    declare -a proteins
    while IFS=',' read -r virus protein; do
      viruses+=("$(echo "$virus" | tr -d '\r')")
      proteins+=("$(echo "$protein" | tr -d '\r')")
    done < <(tail -n +2 "${viral_prots}")
  # Get the total number of entries
    num_entries=${#viruses[@]}

# Directory containing .pdb files
  pdb_dir="pdbs"

# Identify temporary and output directories
  temp="tmp"
    mkdir -p "${temp}"
  output_dir="FATCAT-out"
    mkdir -p "${output_dir}"

# Output .csv file for FATCAT results
  output_file="FATCAT-alignment-Flavi_Prots.csv"
  echo "Protein1,Protein2,ini-len,ini-rmsd,opt-equ,opt-rmsd,chain-rmsd,ali-len,P-value" > "$output_file"


# Nested loops to run FATCAT on virus-protein pairs
  for ((i = 0; i < num_entries; i++)); do
    for ((j = i + 1; j < num_entries; j++)); do
      # Ensure the proteins match
      if [[ "${proteins[$i]}" == "${proteins[$j]}" ]]; then
        # Generate the two variables for the matching proteins
        variable1="${viruses[$i]}-${proteins[$i]}"
        variable2="${viruses[$j]}-${proteins[$j]}"
        echo "Comparing $variable1 and $variable2"
        
        # Set FATCAT output file name
        FATCAT_out="${output_dir}/${variable1}_vs_${variable2}"
        mkdir -p "${FATCAT_out}"
        
        # Set the structure names
        structure1="${pdb_dir}/${variable1}-AF-ranked_0.pdb"
        structure2="${pdb_dir}/${variable2}-AF-ranked_0.pdb"
        
        # Run FATCAT
        FATCAT \
          -p1 "${structure1}" \
          -p2 "${structure2}" \
          -m \
          -ac \
          -t \
          -time \
          -o "${FATCAT_out}/${variable1}_${variable2}"
        
        # Set input file to copy FATCAT .aln file into a .csv file for further analysis
        input_file="${FATCAT_out}/${variable1}_${variable2}.aln"
        
        # Extract values from the FATCAT output
        if [ -f "$input_file" ]; then
          ini_len=$(grep -oP "ini-len \K[0-9]+" "$input_file")
          ini_rmsd=$(grep -oP "ini-rmsd \K[0-9.]+(?= )" "$input_file")
          opt_equ=$(grep -oP "opt-equ \K[0-9]+" "$input_file")
          opt_rmsd=$(grep -oP "opt-rmsd \K[0-9.]+(?= )" "$input_file")
          chain_rmsd=$(grep -oP "chain-rmsd \K[0-9.]+(?= )" "$input_file")
          align_len=$(grep -oP "align-len \K[0-9]+" "$input_file")
          p_value=$(grep -oP "P-value \K[0-9.e+-]+" "$input_file")
          
          # Write extracted data to the .csv output
          echo "$variable1,$variable2,$ini_len,$ini_rmsd,$opt_equ,$opt_rmsd,$chain_rmsd,$align_len,$p_value" >> "$output_file"
        else
          echo "Alignment file not found for $variable1 and $variable2"
        fi
      fi
    done
  done

echo "All done!"
