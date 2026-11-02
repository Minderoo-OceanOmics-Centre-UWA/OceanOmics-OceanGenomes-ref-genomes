module load nextflow/24.10.0
module load singularity/4.1.0-nompi

export NXF_HOME=/scratch/pawsey0964/$USER/.nextflow_home

# Per-run usage reports (trace columns are set in nextflow.config)
NF_TS="$(date +%Y%m%d_%H%M%S)"
NF_INFO_DIR="$PWD/pipeline_info"
mkdir -p "$NF_INFO_DIR"

nextflow run main.nf \
    -profile singularity \
    --input /scratch/pawsey0964/$USER/ref-gen/OceanOmics-OceanGenomes-ref-genomes/assets/samplesheet_20260722.csv \
    --outdir /scratch/pawsey0964/$USER/ref-gen \
    --assembly_mode hifi_hic \
    --scaffolder yahs \
    --buscodb /scratch/pawsey0964/$USER/busco_db/actinopterygii_odb10 \
    --gxdb /scratch/references/Foreign_Contamination_Screening/gxdb \
    --binddir /scratch \
    -c pawsey_profile.config \
    -resume \
    -with-trace "$NF_INFO_DIR/execution_trace_${NF_TS}.txt" \
    -with-report "$NF_INFO_DIR/execution_report_${NF_TS}.html" \
    -with-timeline "$NF_INFO_DIR/execution_timeline_${NF_TS}.html" \
    --tempdir /scratch/pawsey0964/$USER/ref-gen/tmp \
    --bs_config ~/.basespace/default.cfg \
    --sql_config ~/postgresql_details/oceanomics.cfg
