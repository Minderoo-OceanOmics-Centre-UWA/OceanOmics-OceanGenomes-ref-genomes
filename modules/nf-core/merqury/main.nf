process MERQURY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/merqury:1.3--hdfd78af_1':
        'biocontainers/merqury:1.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(meryl_db), path(assembly)
    val(stage)
    val(asmversion)

    output:
    tuple val(meta), path("*.completeness.stats"), emit: stats
    tuple val(meta), path("*.dist_only.hist")    , emit: dist_hist                , optional: true
    tuple val(meta), path("*.spectra-cn.fl.png") , emit: spectra_cn_fl_png        , optional: true
    tuple val(meta), path("*.spectra-cn.hist")   , emit: spectra_cn_hist          , optional: true
    tuple val(meta), path("*.spectra-cn.ln.png") , emit: spectra_cn_ln_png        , optional: true
    tuple val(meta), path("*.spectra-cn.st.png") , emit: spectra_cn_st_png        , optional: true
    tuple val(meta), path("*.spectra-asm.fl.png"), emit: spectra_asm_fl_png       , optional: true
    tuple val(meta), path("*.spectra-asm.hist")  , emit: spectra_asm_hist         , optional: true
    tuple val(meta), path("*.spectra-asm.ln.png"), emit: spectra_asm_ln_png       , optional: true
    tuple val(meta), path("*.spectra-asm.st.png"), emit: spectra_asm_st_png       , optional: true
    tuple val(meta), path("*.qv")                , emit: assembly_qv
    tuple val(meta), path("*.qv")                , emit: scaffold_qv
    tuple val(meta), path("*.hist.ploidy")       , emit: read_ploidy
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 1.3
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi
    # limit meryl to use the assigned number of cores.
    export OMP_NUM_THREADS=$task.cpus

    merqury.sh \\
        $meryl_db \\
        $assembly \\
        ${prefix}.${asmversion}

    # merqury.sh exits 0 even when a step is killed for memory, leaving blank rows.
    # Check every row is filled; exit 137 so the retry with more memory runs.
    n_asm=\$(echo $assembly | wc -w)
    want=\$(( n_asm > 1 ? n_asm + 1 : 1 ))
    got_c=\$(awk -F'\t' '\$5 != ""' ${prefix}.${asmversion}.completeness.stats 2>/dev/null | wc -l)
    got_q=\$(awk -F'\t' '\$4 != ""' ${prefix}.${asmversion}.qv 2>/dev/null | wc -l)
    if [ "\$got_c" -lt "\$want" ] || [ "\$got_q" -lt "\$want" ] || grep -q -i -E 'killed|out of memory|bad_alloc' logs/*.log 2>/dev/null; then
        echo "Merqury output incomplete (completeness \$got_c/\$want, qv \$got_q/\$want rows), probably out of memory" >&2
        exit 137
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merqury: $VERSION
    END_VERSIONS
    """
}
