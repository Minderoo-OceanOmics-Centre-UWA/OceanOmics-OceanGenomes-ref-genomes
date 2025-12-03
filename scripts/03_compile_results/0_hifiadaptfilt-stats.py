import os
import csv

# Define the base directory where the samples are stored
base_dir = "/scratch/pawsey0964/lhuet/ref-gen"  # Change this to the correct path
output_tsv = "hifiadaptorfilt_stats_summary.tsv"

# Dictionary to accumulate totals per sample
sample_totals = {}

# Loop through each sample directory
for sample in os.listdir(base_dir):
    sample_dir = os.path.join(base_dir, sample, '01-data-processing', 'hifiadaptfilt')
    
    if os.path.isdir(sample_dir):
        for root, dirs, files in os.walk(sample_dir):
            for file in files:
                if file.endswith('.stats'):
                    stats_file = os.path.join(root, file)
                    
                    with open(stats_file, 'r') as f:
                        for line in f:
                            if "Number of adapter contaminated ccs reads" in line:
                                contam_reads = line.split(":")[1].split()[0]
                                try:
                                    contam_reads_int = int(contam_reads)
                                except ValueError:
                                    contam_reads_int = 0
                                
                                # Accumulate per sample
                                sample_totals[sample] = sample_totals.get(sample, 0) + contam_reads_int
                                break  # stop after finding the relevant line in this file

# Write output TSV
with open(output_tsv, mode='w', newline='') as tsvfile:
    tsv_writer = csv.writer(tsvfile, delimiter='\t')
    tsv_writer.writerow(['sample', 'contam_reads'])
    
    for sample, total in sorted(sample_totals.items()):
        tsv_writer.writerow([sample, total])

print(f"Summary TSV file generated: {output_tsv}")
