#!/bin/bash

# Path to the file containing the list of genome names
names_list="nom_genomes.txt"

# Check that the file containing the list of genome names exists
if [ ! -f "$names_list" ]; then
    echo "the file $names_list doesn't exist."
    exit 1
fi

mkdir -p .output_fetchMGS

repertory_faa="/PATH/rpob/protein_seq/arc_refseq_cg/2024-02-05_09-01-02/files/modified/"
repertory_fna="/PATH/rpob/cds_genomic_seq/arc_refseq_cg/2024-02-05_09-01-16/files/modified/"
path_general_output="output_fetchMGS/"

# Reading the list of genome names
while IFS= read -r genome; do
    # Building FAA and FNA file paths
    file_faa="${repertory_faa}${genome}_protein_modified.faa"
    file_fna="${repertory_fna}${genome}_cds_from_genomic.fna.gz_modified.fna"

    # Checking the existence of the FNA file
    if [ -e "$fichier_fna" ]; then
        output_folder="${path_general_output}${genome}_output"

        #=fetchMGs.pl
        ./../fetchMGs/fetchMGs.pl -m extraction "$file_faa" -d "$file_fna" -c COG0085 -o "$output_folder"
    else
        echo "Missing .fna file for the genome $genome"
    fi
