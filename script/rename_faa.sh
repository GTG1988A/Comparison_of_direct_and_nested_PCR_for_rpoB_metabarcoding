#!/bin/bash
usage() {
    echo "Usage: $0 <input_directory> <output_directory>"
    echo "  <input_directory>  : Directory containing the .faa.gz files"
    echo "  <output_directory> : Directory where the modified .faa files will be saved"
    exit 1
}

# Checking the number of arguments
if [ "$#" -ne 2 ]; then
    usage
fi

input_dir="$1"
output_dir="$2"

# Check if directories exist
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

# Creation of the output directory if it doesn't exist

if [ ! -d "$output_dir" ]; then
    echo "Output directory '$output_dir' does not exist. Creating it."
    mkdir -p "$output_dir"
fi

#Treatment
for file in "$input_dir"/*.faa.gz; do
    filename=$(basename "$file" .faa.gz)

    zcat "$file" | awk '/^>/ {print $1; next} {print}' > "$output_dir/${filename}_modified.faa"
done

echo "Processing complete. Modified files are saved in '$output_dir'."
