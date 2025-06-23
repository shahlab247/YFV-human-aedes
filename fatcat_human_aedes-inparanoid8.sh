### FATCAT comparison of list of proteins from InParanoid8 after processing the listed output by group
### Reads a .csv file and aligns each protein pair
### Chase Skawinski - 23 JUN 2025
### Copyright (C) 2025 Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu

#!/bin/bash

## Identify Human-Aedes paralogs file
  paralog_file="inParanoid8-Human_Aedes-Homologs.csv"
    # Add a newline to the end of the CSV file if it is missing
    if [ -s "${paralog_file}" ] && [ "$(tail -c1 "${paralog_file}")" != "" ]; then
      echo "" >> "${paralog_file}"
    fi

## Set paths to proteome .tar files
  aedes_tar="PATH_TO_TAR"
  human_tar="PATH_TO_TAR"

## Generate .csv file for alignment outputs
  output_file="FATCAT-alignment-Human_Aedes.csv"
    # Add the relevant headers to the .csv file
    echo "Group,Human Protein,Aedes Protein,ini-len,ini-rmsd,opt-equ,opt-rmsd,chain-rmsd,ali-len,P-value" > $output_file

# Set path to a temp directory
  tmp="tmp"
  mkdir -p "${tmp}"

# Initialize variables
  current_group=""
  aedes_proteins=()
  homo_proteins=()


## Function to extract .pdb files from a .tar archive if they exist
extract_pdb_files() {
  local tar_file=$1
  local protein_id=$2
  local temp_dir=$3
  
  # Ensure the temp directory exists
    mkdir -p "${temp_dir}"

  # Attempt extraction of .pdb.gz and .pdb files
    echo "  Searching for files..."
    tar -xf "${tar_file}" -C "${temp_dir}" --wildcards --no-anchored *"${protein_id}"*.pdb*

  # Check if .pdb.gz files are present and decompress them
    if ls "${temp_dir}"/*"${protein_id}"*.pdb.gz >/dev/null 2>&1; then
      echo "  Decompressing .pdb.gz files in $temp_dir"
      gunzip "$temp_dir"/*"${protein_id}"*.pdb.gz
      if [ $? -ne 0 ]; then
        echo "  Error decompressing .pdb.gz files"
        return 1
      fi
    fi

  # Verify and list extracted .pdb files
    extracted_files=($(ls "${temp_dir}"/*"${protein_id}"*.pdb 2>/dev/null))
    if [ ${#extracted_files[@]} -eq 0 ]; then
      echo "  No .pdb file found for $protein_id in ${tar_file}"
      return 1
    else
      echo "  Extracted files: ${extracted_files[@]}"
      return 0
    fi
}


## Read the .csv file and extract the .pdb file(s) for each homolog
while IFS=, read -r group bitscore organism inparalog_score protein_id bootstrap_support protein_name protein_length; do
  if [ "${organism}" == "Aedes aegypti" ]; then
    echo "Beginning extraction of '${protein_id}' from '${aedes_tar}'"
    extract_pdb_files "${aedes_tar}" "${protein_id}" "${tmp}"
    echo "Extraction of '${protein_id}' from '${aedes_tar}' complete"
  elif [ "${organism}" == "Homo sapiens" ]; then
    echo "Beginning extraction of ${protein_id} from ${human_tar}"
    extract_pdb_files "${human_tar}" "${protein_id}" "${tmp}"
    echo "Extraction of '${protein_id}' from '${human_tar}' complete"
  else
    echo "Organism ${organism} not found."
  fi
done < <(tail -n +2 "${paralog_file}") # Iterates to the next line for the while function


## Function to run FATCAT on corresponding protein pairs
run_fatcat_for_group() {
  if [ -n "${current_group}" ]; then
    echo "  Running FATCAT for group: ${current_group}"
    for protein_aedes in "${aedes_proteins[@]}"; do
      for protein_homo in "${homo_proteins[@]}"; do
      # Search for all .pdb files related to each protein ID
        echo "    Aedes aegypti Protein ID: $protein_aedes"
        echo "    Homo sapiens Protein ID: $protein_homo"
        aedes_pdb_files=($(ls "${tmp}"/*"${protein_aedes}"*.pdb 2>/dev/null))
        homo_pdb_files=($(ls "${tmp}"/*"${protein_homo}"*.pdb 2>/dev/null))
        
      # Ensure the .pdb files are present to process and print to output if they aren't
        if [ ${#homo_pdb_files[@]} -eq 0 ] && [ ${#aedes_pdb_files[@]} -eq 0 ]; then
          echo "    No .pdb files found for '$protein_homo' or '$protein_aedes'"
          echo "$current_group,MISSING-$protein_homo,MISSING-$protein_aedes" >> "${output_file}"
        elif [ ${#homo_pdb_files[@]} -eq 0 ]; then
          echo "    No .pdb file found for '$protein_homo'"
          echo "$current_group,MISSING-$protein_homo,$protein_aedes" >> "${output_file}"
        elif [ ${#aedes_pdb_files[@]} -eq 0 ]; then
          echo "    No .pdb file found for '$protein_aedes'"
          echo "$current_group,$protein_homo,MISSING-$protein_aedes" >> "${output_file}"
          continue
        fi
        
      # Iterate over each .pdb file combination
        for aedes_pdb in "${aedes_pdb_files[@]}"; do
          for homo_pdb in "${homo_pdb_files[@]}"; do
            
          #Extract the basename of each .pdb file (without extension)
            aedes_base=$(basename "${aedes_pdb}" .pdb)
            homo_base=$(basename "${homo_pdb}" .pdb)
            echo "    Human base file name is $homo_base"
            echo "    Aedes base file name is $aedes_base" 

          # Create output directory for FATCAT
            FATCAT_out="FATCAT-out/${current_group}_${homo_base}_${aedes_base}"
            mkdir -p "${FATCAT_out}"

          # Run FATCAT
            echo "Running FATCAT with ${homo_base} and ${aedes_base}"
            FATCAT \
              -p1 "${tmp}"/"${homo_base}".pdb \
              -p2 "${tmp}"/"${aedes_base}".pdb \
              -m \
              -ac \
              -t \
              -time \
              -o "${FATCAT_out}"/"${homo_base}"_"${aedes_base}"

          # Check if FATCAT output exists
            if [ -f "${FATCAT_out}/${homo_base}_${aedes_base}.aln" ]; then
              input_file="${FATCAT_out}/${homo_base}_${aedes_base}.aln"

            # Extract values from the FATCAT output
              ini_len=$(grep -oP "ini-len \K[0-9]+" "${input_file}")
              ini_rmsd=$(grep -oP "ini-rmsd \K[0-9.]+(?= )" "${input_file}")
              opt_equ=$(grep -oP "opt-equ \K[0-9]+" "${input_file}")
              opt_rmsd=$(grep -oP "opt-rmsd \K[0-9.]+(?= )" "${input_file}")
              chain_rmsd=$(grep -oP "chain-rmsd \K[0-9.]+(?= )" "${input_file}")
              align_len=$(grep -oP "align-len \K[0-9]+" "${input_file}")
              p_value=$(grep -oP "P-value \K[0-9.e+-]+" "${input_file}")

            # Write the extracted data to the CSV file
              echo "$current_group,$homo_base,$aedes_base,$ini_len,$ini_rmsd,$opt_equ,$opt_rmsd,$chain_rmsd,$align_len,$p_value" >> "${output_file}"
              echo "  Alignment .csv file '$output_file' updated with data from $input_file."
            else
              echo "  FATCAT output for ${protein_homo}_${protein_aedes} not found."
            fi
          done
        done
      done
    done  
  fi
}


## Read the homolog .csv file and process each group
while IFS=, read -r group bitscore organism inparalog_score protein_id bootstrap_support protein_name protein_length; do
  [ -z "${group}" ] && continue
  echo "Processing line: group=$group, organism=$organism, protein_id=$protein_id"
  
  # Check if in a new group
  if [ "${current_group}" != "${group}" ] && [ -n "${current_group}" ]; then
  
  # Run FATCAT for the previous group
    echo "About to process group: $current_group"
    run_fatcat_for_group
  # Reset the variables for the new group
    aedes_proteins=()
    homo_proteins=()
  fi
  
  # Update the new group
    current_group="${group}"

  # Collect protein IDs based on the organism
  if [ "${organism}" == "Aedes aegypti" ]; then
    aedes_proteins+=("${protein_id}")
  elif [ "${organism}" == "Homo sapiens" ]; then
    homo_proteins+=("${protein_id}")
  fi
  
done < <(tail -n +2 "${paralog_file}") # Iterates to the next line for the while function


## Ensure FATCAT runs for the final group after the loop finishes
if [ -n "${current_group}" ]; then
  echo "Processing final group: ${current_group}"
  run_fatcat_for_group
fi

## Delete the extracted .pdb files to clear up space once we're all done!
rm -f "${tmp}"/*.pdb*
