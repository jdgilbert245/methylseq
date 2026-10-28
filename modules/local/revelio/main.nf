process REVELIO {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py310h4a09ff2_1'
        : 'quay.io/biocontainers/pysam:0.24.0--py310h4a09ff2_1'}"

    input:
    tuple val(meta), path(bam_calmd), path(bai_calmd)

    output:
    tuple val(meta), path("*.masked.bam"), emit: masked_bam

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    revelio.py \\
        -Q \\
        -T ${task.cpus} \\
        -t . \\
        ${args} \\
        ${bam_calmd} \\
        ${prefix}.masked.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.masked.bam
    """
}
