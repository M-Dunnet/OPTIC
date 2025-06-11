#!/bin/bash
"""
# This script moves files listed in a TSV file from an input directory to an output directory.
# Its purpose is to organise files seperated by the hierarchical clustering function of OPTIC
# Usage: ./move_files.sh <tsv_file> <input_directory> <output_directory>
"""


if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <tsv_file> <input_directory> <output_directory>"
    exit 1
fi

tsv_file="$1"
input_dir="$2"
output_dir="$3"

# Get the number of columns (headers count)
num_columns=$(head -n 1 "$tsv_file" | awk -F'\t' '{print NF}')

# Process each column separately
for (( col=1; col<=num_columns; col++ )); do
    # Extract the column, skipping the header, and process it
    awk -v col="$col" 'BEGIN {FS=OFS="\t"} NR > 1 {print $col}' "$tsv_file" | while IFS= read -r file; do
        if [ -n "$file" ]; then
            file_path="${input_dir}/${file}"
            if [ -e "$file_path" ]; then
                mkdir -p "${output_dir}/Column_${col}_samples/"
                cp "$file_path"* "${output_dir}/Column_${col}_samples/"
            else
                echo "Warning: File $file_path not found"
            fi
        fi
    done
done
