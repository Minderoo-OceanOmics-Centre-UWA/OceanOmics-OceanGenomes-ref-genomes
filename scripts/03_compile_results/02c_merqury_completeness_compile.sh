#!/bin/bash

output_file="merqury.completeness.stats.tsv"
base_dir="/scratch/pawsey0964/$USER/ref-gen/OG*"

# Write header
echo -e "sample\tsolid_k_mers\ttotal_k_mers\tcompleteness" > "$output_file"

# Find all .hifiasm.completeness.stats files
completeness_files=$(find $base_dir -name "*.hifiasm.completeness.stats")

for file in $completeness_files; do
    sample_base=$(basename "$file" .hifiasm.completeness.stats)

    while IFS=$'\t' read -r sample _ solid_kmers total_kmers completeness; do
        # Handle dual-haplotype entry
        if [[ "$sample" == "both" || "$sample" == "Both" ]]; then
            sample="${sample_base}.hifiasm.dual"
        fi

        echo -e "${sample}\t${solid_kmers}\t${total_kmers}\t${completeness}" >> "$output_file"
    done < "$file"
done