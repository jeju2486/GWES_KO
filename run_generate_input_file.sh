#!/bin/bash

# Define the input file and destination directory
input_file="input_100.txt"  # Change this to the path of your input file containing .fas paths
output_file="input_100_gff.txt"
destination_dir="/data/biol-micro-genomics/kell7366/minimap/ref_fas/ref_aureus/100_gffs"  # Change this to your desired destination directory

# Create the destination directory if it doesn't exist
mkdir -p "$destination_dir"

# Clear the output file if it exists
> "$output_file"

# Loop through each line in the input file
while IFS= read -r line; do
  # Replace 'negative' with 'negative_gff' and 'positive' with 'positive_gff' and change .fas to _prokka.gff
  modified_line=$(echo "$line" | sed 's/\/negative\//\/negative_gff\//; s/\/positive\//\/positive_gff\//; s/\.fas/_prokka.gff/')
  # Append the modified line to the output file
  echo "$modified_line" >> "$output_file"
  
  # Copy the file to the destination directory if it exists
  if [[ -f "$modified_line" ]]; then
    cp "$modified_line" "$destination_dir/"
    echo "Copied $modified_line to $destination_dir/"
  else
    echo "File $modified_line does not exist."
  fi
done < "$input_file"

echo "Process completed. Check $output_file for the modified paths and $destination_dir for the copied files."
