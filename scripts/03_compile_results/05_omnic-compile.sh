#!/bin/bash
# Compile OMNIC mapping stats with stage information
base_dir="/scratch/pawsey0964/$USER/ref-gen" 
# Print header
echo -e "sample\tstage\thaplotype\ttotal\ttotal_unmapped\ttotal_single_sided_mapped\ttotal_mapped\ttotal_dups\ttotal_nodups\tcis\ttrans" > final_omnic_stats_report.txt

# Find all *.stats.txt files
stats_files=$(find $base_dir/OG* -name "*.stats.txt")

# Process each stats.txt file
for file in $stats_files; do

    # Extract sample name (OG*) and haplotype (hap1, hap2, dual)
    filename=$(basename "$file")
    sample=$(echo "$filename" | cut -d'.' -f1)
    
    # Extract the "stage" — the parent directory two levels up (04-scaffolding or 05-decontamination)
    stage=$(basename "$(dirname "$(dirname "$file")")")
    
    if [[ "$filename" == *hap1* ]]; then
        haplotype="hap1"
    elif [[ "$filename" == *hap2* ]]; then
        haplotype="hap2"
    elif [[ "$filename" == *dual* ]]; then
        haplotype="dual"
    else
        haplotype="unknown"
    fi

    # Extract the numbers
    total=$(awk '/^total[[:space:]]/ {print $2}' "$file")
    total_unmapped=$(awk '/^total_unmapped/ {print $2}' "$file")
    total_single_sided_mapped=$(awk '/^total_single_sided_mapped/ {print $2}' "$file")
    total_mapped=$(awk '/^total_mapped/ {print $2}' "$file")
    total_dups=$(awk '/^total_dups/ {print $2}' "$file")
    total_nodups=$(awk '/^total_nodups/ {print $2}' "$file")
    cis=$(awk '/^cis[[:space:]]/ {print $2}' "$file")
    trans=$(awk '/^trans[[:space:]]/ {print $2}' "$file")

    # Print to output
    echo -e "${sample}\t${stage}\t${haplotype}\t${total}\t${total_unmapped}\t${total_single_sided_mapped}\t${total_mapped}\t${total_dups}\t${total_nodups}\t${cis}\t${trans}" >> final_omnic_stats_report.txt

done
