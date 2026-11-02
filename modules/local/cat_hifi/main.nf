process CAT_HIFI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(files)

    output:
    tuple val(meta), path("*.hifi.cat.fastq.gz"), emit: hifi_cat
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

        cat \\
            $args \\
            *fastq.gz \\
            > ${prefix}.hifi.cat.fastq.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busybox: \$( busybox | sed '1!d; s/BusyBox v//; s/ .*//' )
    END_VERSIONS
    """
}
