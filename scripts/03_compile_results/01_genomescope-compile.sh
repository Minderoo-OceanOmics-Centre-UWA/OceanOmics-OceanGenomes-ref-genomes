#!/bin/bash

base_dir="/scratch/pawsey0964/$USER/ref-gen" 
# Define the output file and create the column headings
TSV="genomescope_compiled_results.tsv"
echo sample,homozygosity,heterozygosity,genomesize,repeatsize,uniquesize,modelfit,errorrate | sed 's/,/\t/g' | tee $TSV


for GSCOPE in $base_dir/OG*/02-kmer-profiling/genomescope2/*_genomescope_summary.txt; do
PREFIX=$(basename $GSCOPE | awk -F "-genomescope" '{print $1;}')
sed -E 's/ [ ]+/\t/g' $GSCOPE | awk -F '\t' '{print $3;}' | tail -n 8 | awk '{print $1;}' | tr '\n' '\t' | sed 's/\t$/\n/' | sed "s/max/$PREFIX/" | tee -a $TSV
done
column -t $TSV
