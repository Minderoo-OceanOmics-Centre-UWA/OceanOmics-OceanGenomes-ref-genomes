process OMNIC {
    tag "$meta.id"
    label 'process_medium'

    //conda "${moduleDir}/environment.yml"
    container "docker://sawtooth01/omnic:v0.01"

    input:
    tuple val(meta), path(reads), path(assembly)
    val(haplotype)
    val(tempdir)

    output:
    tuple val(meta), path("*.stats.txt"), emit: omnic_stats
    tuple val(meta), path("*.mapped.PT.bam"), emit: omnic_bam
    tuple val(meta), path("*.mapped.PT.bam.csi"), emit: omnic_bai
    tuple val(meta), path("*fai"), emit: omnic_fai
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    export PATH=\$PATH:/opt/conda/envs/pairtools/bin

    samtools faidx ${assembly}
    cut -f1,2 ${assembly}.fai > "${meta.id}.${haplotype}.genome"
    bwa index ${assembly}
    bwa mem -5SP -T0 -t24 ${assembly} ${reads} \\
    | pairtools parse \\
        --min-mapq 0 \\
        --walks-policy all \\
        --nproc-in 24 \\
        --nproc-out 24 \\
        --chroms-path "${meta.id}.${haplotype}.genome" \\
    | pairtools sort \\
        --nproc 24 \\
        --tmpdir ${tempdir} \\
    | pairtools dedup \\
        --nproc-in 24 \\
        --nproc-out 24 \\
        --mark-dups \\
        --output-stats "${meta.id}.${haplotype}.stats.txt" \\
    | pairtools split \\
        --nproc-in 24 \\
        --nproc-out 24 \\
        --output-pairs /dev/null \\
        --output-sam - \\
    | samtools sort -@32 -T "${tempdir}/${meta.id}_temp.bam" -o "${meta.id}.${haplotype}.mapped.PT.bam" -

    samtools index --csi "${meta.id}.${haplotype}.mapped.PT.bam"

    rm -f "${meta.id}.${haplotype}.genome" \\
          ${assembly}.amb ${assembly}.ann ${assembly}.bwt ${assembly}.pac ${assembly}.sa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed '/^Version:/!d; s/Version: //')
        pairtools: \$(pairtools --version 2>&1 | sed 's/.*version //')
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch "${meta.id}.stats.txt"
    touch "${meta.id}.mapped.PT.bam"
    touch "${meta.id}.mapped.PT.bam.csi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | sed '/^Version:/!d; s/Version: //')
        pairtools: \$(pairtools --version 2>&1 | sed 's/.*version //')
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
