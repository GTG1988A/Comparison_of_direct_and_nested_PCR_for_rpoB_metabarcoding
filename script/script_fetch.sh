#!/bin/bash

# Function to display help message
usage() {
    echo "Usage: $0 genome_names_file global_output_path fna_directory faa_directory"
    echo "  genome_names_file   : File containing the list of genome names"
    echo "  global_output_path  : Directory where the output will be saved"
    echo "  fna_directory       : Directory containing the modified .fna files"
    echo "  faa_directory       : Directory containing the modified .faa files"
    exit 1
}

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    usage
fi

genome_names_file=$1
global_output_path=$2
fna_directory=$3
faa_directory=$4

# Check if the genome names file exists
if [ ! -f "$genome_names_file" ]; then
    echo "The file $genome_names_file does not exist."
    exit 1
fi

# Create the global output directory if it doesn't exist
mkdir -p "$global_output_path"

# Check if the FNA and FAA directories exist
if [ ! -d "$fna_directory" ]; then
    echo "The FNA directory $fna_directory does not exist."
    exit 1
fi

if [ ! -d "$faa_directory" ]; then
    echo "The FAA directory $faa_directory does not exist."
    exit 1
fi

# Read the genome names file
while IFS= read -r genome; do
    # Construct the paths for the FAA and FNA files
    faa_file="${faa_directory}/${genome}_protein_modified.faa"
    fna_file="${fna_directory}/${genome}_cds_from_genomic_modified.fna"

    # Check if the FNA file exists
    if [ -e "$fna_file" ]; then
        # Create the output directory
        output_directory="${global_output_path}/${genome}_output"

        # Call fetchMGs.pl with the paths of the FAA and FNA files
        ./fetchMGs/fetchMGs.pl -m extraction "$faa_file" -d "$fna_file" -c COG0085 -o "$output_directory"
    else
        echo "Missing .fna file for genome $genome"
    fi
done < "$genome_names_file"
