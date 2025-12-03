module load nextflow/24.10.0
module load singularity/4.1.0-nompi

nextflow run main.nf \
    -profile singularity \
    --input /scratch/pawsey0964/lhuet/ref-gen/OceanOmics-OceanGenomes-ref-genomes/assets/samplesheet.csv \
    --outdir /scratch/pawsey0964/lhuet/ref-gen \
    --scaffolder yahs \
    --buscodb /scratch/references/busco_db/actinopterygii_odb10 \
    --gxdb /scratch/references/Foreign_Contamination_Screening/gxdb \
    --binddir /scratch \
    -c pawsey_profile.config \
    -resume \
    --tempdir /scratch/pawsey0964/lhuet/ref-gen/tmp \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg
