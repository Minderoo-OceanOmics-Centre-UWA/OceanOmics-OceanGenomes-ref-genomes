process PRETEXTGRAPH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pretextgraph:0.0.7--h9948957_1':
        'biocontainers/pretextgraph:0.0.7--h9948957_1' }"

    input:
    tuple val(meta), path(pretext_map), path(bedgraph)
    val(label)

    output:

    tuple val(meta), path(pretext_map)  , emit: pretext_graph
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $bedgraph | \\
    PretextGraph \\
        -i $pretext_map \\
        -n $label \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextgraph: \$(PretextGraph | grep "Version" | sed 's/.*Version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pretextgraph: \$(PretextGraph | grep "Version" | sed 's/.*Version //')
    END_VERSIONS
    """
}
