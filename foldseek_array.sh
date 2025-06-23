### Foldseek comparison of list of proteins to a database - either human or Aedes aegypti
### Reads a .csv file and aligns each protein against the other organism's database
### Change the paths and directories for correct organism based on comparison needed
### Chase Skawinski - 23 JUN 2025
### Copyright (C) 2025 Shah Lab, UC Davis
### GNU General Public License v3.0
### email: prsshah@ucdavis.edu


#!/bin/bash

# Foldseek database paths
  afaedes="PATH_TO_DIRECTORY"
  afhuman="PATH_TO_DIRECTORY"

# Protein structure .tar paths
  aedes_tar="PATH_TO_DIRECTORY"
  human_tar="PATH_TO_DIRECTORY"
  
# Set which Foldseek database and .tar file you want to use
  fold_db="${afhuman}"
  tar_file="${aedes_tar}"

# Protein list path
  protein_file="aedes_proteins_list.txt"
    # Add a newline to the end of the file if it's missing
    if [ -s "${protein_file}" ] && [ "$(tail -c1 "${protein_file}")" != "" ]; then
      echo "" >> "${protein_file}"
    fi
  
# Output path
  output_dir="Aedes_to_Human_full"
    mkdir -p "${output_dir}"

# Make a temp directory
  temp_dir="temp"
    mkdir -p "${temp_dir}"

# Create an output .csv file
  output_file="foldseek_results.csv"
  
  # Write the headers to the .csv file
  output_headers="query,qlen,qcov,target,tlen,tcov,alnlen,bits,mismatch,evalue,alntmscore,qtmscore,ttmscore,lddt,prob"
  echo "${output_headers}" > "${output_file}"


## Run Foldseek iteratively with proteins in the list
while IFS= read -r protein_file; do
  protein_id=${protein_file} 
  echo "Processing ${protein_id}..."
  
  # Extract the .pdb file from the .tar
  tar -xf "${aedes_tar}" -C "${temp_dir}" --wildcards --no-anchored *"${protein_id}"*.pdb*
  
  # Check if a .pdb.gz file is present and decompress it
  if ls "${temp_dir}"/*"${protein_id}"*.pdb.gz >/dev/null 2>&1; then
    echo "  Decompressing .pdb.gz files in ${temp_dir}"
    gunzip "${temp_dir}"/*"${protein_id}"*.pdb.gz
    if [ $? -ne 0 ]; then
      echo "  Error decompressing .pdb.gz files"
      continue
    fi
  fi

  # Verify and list the extracted .pdb file(s)
  extracted_files=($(ls "${temp_dir}"/*"${protein_id}"*.pdb 2>/dev/null))
  if [ ${#extracted_files[@]} -eq 0 ]; then
    echo "    No .pdb file found for ${protein_id} in ${tar_file}"
    continue
  else
    echo "    Extracted files: ${extracted_files[@]}"
  fi

  # Use a for-loop to iterate over protein comparisons if multiple .pdb files are present 
  pdb_files=($(ls "${temp_dir}"/*"${protein_id}"*.pdb 2>/dev/null))
  for pdb in "${pdb_files[@]}"; do
    echo "Running Foldseek on ${pdb}"
      
    # Extract the basename of each file
      protein_base=$(basename "${pdb}" .pdb)
      echo "  Base file name is ${protein_base}"
      query_output="${output_dir}/${protein_base}_ali.txt"
         
    # Run Foldseek
      foldseek \
        easy-search \
        "${pdb}" \
        "${fold_db}" \
        "${query_output}" \
        tmpFolder \
        --format-mode 0 \
        --format-output "${output_headers}" \
        -s 10 \
        -e 1E-10
             
    if [ $? -ne 0 ]; then
      echo "  Foldseek failed on ${pdb}"
      continue
    fi

    # Add results to the output .csv file
    tail -n +2 "${query_output}" | tr '\t' ',' >> "${output_file}"

    echo "  Foldseek completed for ${pdb}"
  done

    # Clean up temporary files
      rm -rf "${temp_dir}"/*
done < "${protein_file}"

echo "All Foldseek searches complete."
