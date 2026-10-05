process HIFIADAPTERFILT {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/biocontainers/hifiadapterfilt:2.0.0_cv2'

    input:
    tuple val(meta), path(reads_dir)

    output:
    tuple val(meta), path("$meta.id/*.fastq.gz"), emit: reads
    tuple val(meta), path("$meta.id/*.stats")   , emit: stats
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $meta.id

    if ls $reads_dir/*.bam 1> /dev/null 2>&1; then
        ln -s $reads_dir/*bam .

        bash pbadapterfilt.sh \\
            $args \\
            -t $task.cpus \\
            -o $meta.id
    else
        cp $reads_dir/*fastq.gz $meta.id

        touch $meta.id/dummy.stats
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HiFiAdapterFilt: 2.0.0
        blastn: \$(blastn -version | sed '1!d; s/blastn: //; s/+\$//')
        bamtools: \$(bamtools --version | sed '/^bamtools /!d; s/bamtools //')
    END_VERSIONS
    """
}
