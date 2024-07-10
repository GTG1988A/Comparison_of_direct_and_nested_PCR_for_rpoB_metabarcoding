#!/bin/bash

# Fonction pour afficher l'aide
usage() {
    echo "Usage: $0 <input_directory> <output_directory>"
    echo "  <input_directory>  : Directory containing the .fna.gz files"
    echo "  <output_directory> : Directory where the modified .fna files will be saved"
    exit 1
}

# Vérification du nombre d'arguments
if [ "$#" -ne 2 ]; then
    usage
fi

# Répertoires d'entrée et de sortie
input_dir="$1"
output_dir="$2"

# Vérification si le répertoire d'entrée existe
if [ ! -d "$input_dir" ]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

# Création du répertoire de sortie s'il n'existe pas
if [ ! -d "$output_dir" ]; then
    echo "Output directory '$output_dir' does not exist. Creating it."
    mkdir -p "$output_dir"
fi

# Traitement des fichiers
for file in "$input_dir"/*.fna.gz; do
    filename=$(basename "$file" .fna.gz)

    zcat "$file" | awk '/^>/ {
        if (match($0, /\[protein_id=([^]]*)\]/, arr)) {
            print ">" arr[1];
            skip_seq = 0;
        } else {
            skip_seq = 1;
        }
        next
    } !skip_seq {print}' > "$output_dir/${filename}_modified.fna"
done

echo "Processing complete. Modified files are saved in '$output_dir'."
