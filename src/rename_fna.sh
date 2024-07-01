#!/bin/bash

for file in PATH/fna/arc_refseq_cg/2024-02-29_13-59-50/files/*.fna.gz; do
    filename=$(basename "$file" .fna.gz)

    zcat "$file" | awk '/^>/ {
        if (match($0, /\[protein_id=([^]]*)\]/, arr)) {
            print ">" arr[1];
            skip_seq = 0;
        } else {
            skip_seq = 1;
        }
        next
    } !skip_seq {print}' > "PATH/fna/arc_refseq_cg/2024-02-29_13-59-50/files/modified/${filename}_modified.fna"
done
