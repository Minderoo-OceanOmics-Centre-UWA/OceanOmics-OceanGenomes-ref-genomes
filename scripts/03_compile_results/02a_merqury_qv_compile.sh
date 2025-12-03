output_file="merqury.qv.stats.tsv"
base_dir="/scratch/pawsey0964/$USER/ref-gen/OG*" 
echo -e "sample\tunique_k_mers_assembly\tk_mers_total\tqv\terror" > "$output_file"

# Find all .hifiasm.qv files
completeness_files=$(find $base_dir -name "*.hifiasm.qv")

for file in $completeness_files; do
    sample_base=$(basename "$file" .hifiasm.qv)
    
    while IFS=$'\t' read -r sample unique_k_mers k_mers_total qv error; do
        # Trim whitespace
        trimmed_sample=$(echo "$sample" | xargs)

        if [[ "$trimmed_sample" == "Both" ]]; then
            new_sample="${sample_base}.hifiasm.dual"
            echo -e "${new_sample}\t${unique_k_mers}\t${k_mers_total}\t${qv}\t${error}" >> "$output_file"
        else
            echo -e "${sample}\t${unique_k_mers}\t${k_mers_total}\t${qv}\t${error}" >> "$output_file"
        fi
    done < "$file"
done
