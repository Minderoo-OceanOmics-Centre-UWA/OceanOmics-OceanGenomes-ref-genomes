process TIARA_TIARA {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tiara:1.0.3' :
        'biocontainers/tiara:1.0.3' }"

    input:
    tuple val(meta), path(fasta)
    val(haplotype)

    output:
    tuple val(meta), path("${prefix}*.{txt,txt.gz}"), emit: classifications
    tuple val(meta), path("${removal}")       , emit: contig_removal
    tuple val(meta), path("${review}")        , emit: kept_for_review
    tuple val(meta), path("log_*.{txt,txt.gz}")     , emit: log
    tuple val(meta), path("*.{fasta,fasta.gz}")     , emit: fasta, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tiara_report="${prefix}${haplotype}.txt"
    summary="${prefix}${haplotype}.tiara_filter_summary.txt"
    removal="${prefix}${haplotype}.tiara.contig_removal.txt"
    review="${prefix}${haplotype}.tiara.kept_for_review.tsv"
    def keep_min = task.ext.prok_keep_min_len != null ? task.ext.prok_keep_min_len : 0
    def cat_cmd = fasta.name.endsWith('.gz') ? 'gzip -dc' : 'cat'
    """
    tiara -i ${fasta} \
        -o ${prefix}${haplotype}.txt \
        --threads ${task.cpus} \
        ${args}

    find . -name "*_${fasta}*" -exec sh -c 'file=`basename {}`; mv "\$file" "\${file%%_*}_${prefix}_${haplotype}.fasta"' \\;

    # Sequence lengths (the Tiara report has none, which made the old bp column always 0)
    ${cat_cmd} ${fasta} | awk '/^>/ {if (n != "") print n"\\t"l; n = substr(\$1, 2); l = 0; next} {l += length(\$0)} END {if (n != "") print n"\\t"l}' > seq_lengths.tsv

    # One row per non-eukaryotic call: id, class, length, action.
    # Class is the second-stage call for organelles, the first-stage call otherwise.
    # Mitochondrion/plastid calls are always removed. Prokaryote calls (prokarya/bacteria/archaea) of
    # at least ${keep_min} bp are kept and listed for review: these sequences have already passed
    # FCS-GX, and Tiara calls large repeat-rich fish scaffolds "archaea" (OG866 hic2: 8 scaffolds,
    # 0.16-1.48 Mb, at normal HiFi depth). keep_min = 0 restores the old remove-everything behaviour.
    awk -F'\\t' -v min=${keep_min} 'BEGIN {OFS = "\\t"}
        NR == FNR {len[\$1] = \$2; next}
        FNR == 1  {next}
        {
            c = (\$2 == "organelle") ? \$3 : \$2
            L = len[\$1] + 0
            if (c == "mitochondrion" || c == "plastid") print \$1, c, L, "remove"
            else if (c == "prokarya" || c == "bacteria" || c == "archaea") print \$1, c, L, ((min > 0 && L >= min) ? "keep_review" : "remove")
        }' seq_lengths.tsv "$tiara_report" > tiara_calls.tsv

    awk -F'\\t' '\$4 == "remove" {print \$1}' tiara_calls.tsv > "$removal"

    printf "sequence_id\\ttiara_class\\tlength_bp\\treason\\n" > "$review"
    awk -F'\\t' -v min=${keep_min} 'BEGIN {OFS = "\\t"} \$4 == "keep_review" {print \$1, \$2, \$3, "prokaryote call >= " min " bp after FCS-GX; kept, check HiFi depth/BUSCO before removing"}' tiara_calls.tsv >> "$review"

    # Summary: all calls per class (num_contigs, bp), then what was removed and what was kept for review
    printf "Category\\tnum_contigs\\tbp\\n" > "$summary"
    for c in mitochondrion plastid prokarya bacteria archaea; do
        awk -F'\\t' -v c=\$c 'BEGIN {n = 0; b = 0} \$2 == c {n++; b += \$3} END {printf "%s\\t%d\\t%d\\n", toupper(substr(c, 1, 1)) substr(c, 2), n, b}' tiara_calls.tsv >> "$summary"
    done
    awk -F'\\t' 'BEGIN {n = 0; b = 0} \$4 == "remove"      {n++; b += \$3} END {printf "Removed\\t%d\\t%d\\n", n, b}' tiara_calls.tsv >> "$summary"
    awk -F'\\t' 'BEGIN {n = 0; b = 0} \$4 == "keep_review" {n++; b += \$3} END {printf "Kept_for_review\\t%d\\t%d\\n", n, b}' tiara_calls.tsv >> "$summary"
    rm -f seq_lengths.tsv tiara_calls.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiara: ${VERSION}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    removal="${prefix}${haplotype}.tiara.contig_removal.txt"
    review="${prefix}${haplotype}.tiara.kept_for_review.tsv"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.out.txt
    touch log_${prefix}.out.txt
    touch bacteria_${prefix}.fasta
    touch ${prefix}${haplotype}.tiara.contig_removal.txt
    touch ${prefix}${haplotype}.tiara.kept_for_review.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tiara: ${VERSION}
    END_VERSIONS
    """
}
