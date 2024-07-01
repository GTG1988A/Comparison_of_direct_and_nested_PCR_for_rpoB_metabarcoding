#!/bin/bash

for file in /path_to_data/protein_seq/arc_refseq_cg/2024-02-05_09-01-02/files/*.faa.gz; do

    filename=$(basename "$file" .faa.gz)

    zcat "$file" | awk '/^>/ {print $1; next} {print}' > "/path_to_data/protein_seq/arc_refseq_cg/2024-02-05_09-01-02/files/modified/${filename}_modified.faa"
done
