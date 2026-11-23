import os
import sys
import csv

# Accept base directory as CLI arg; fall back to standard scratch layout for the current user
base_dir = sys.argv[1] if len(sys.argv) > 1 else f"/scratch/pawsey0964/{os.environ['USER']}/ref-gen"
output_tsv = "hifiadaptorfilt_stats_summary.tsv"
# One row per .stats file (one per cell BAM): 0a pushes these per sequencing run
# (ref_hifi_adapter_filter_run, genomes_db_review migration 093).
files_tsv = "hifiadaptorfilt_stats_files.tsv"

# Dictionary to accumulate totals per sample
sample_totals = {}
stats_files = []

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
                                stats_files.append((sample, os.path.relpath(stats_file, base_dir), contam_reads_int))
                                break  # stop after finding the relevant line in this file

# Write output TSV
with open(output_tsv, mode='w', newline='') as tsvfile:
    tsv_writer = csv.writer(tsvfile, delimiter='\t')
    tsv_writer.writerow(['sample', 'contam_reads'])
    
    for sample, total in sorted(sample_totals.items()):
        tsv_writer.writerow([sample, total])

with open(files_tsv, mode='w', newline='') as tsvfile:
    tsv_writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
    tsv_writer.writerow(['sample', 'path', 'contam_reads'])
    for row in sorted(stats_files):
        tsv_writer.writerow(row)

print(f"Summary TSV file generated: {output_tsv}")
print(f"Per-file TSV generated: {files_tsv}")
