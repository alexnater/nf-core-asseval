process PANDEPTH {
    tag "$meta.id"
    label 'process_medium'

    // FIXME Conda is not supported at the moment
    container "${projectDir}/assets/containers/pandepth_2.25.sif"

    input:
    tuple val(meta) , path(bam), path(bai), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(gff)

    output:
    tuple val(meta), path("*.stat.gz"), emit: stats
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PANDEPTH module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "-r ${fasta}" : ""
    def annotation = gff ? "-g ${gff}" : ""
    def interval = bed ? "-b ${bed}" : ""
    def VERSION = '2.25'

    """
    pandepth \\
        -t $task.cpus \\
        $interval \\
        $reference \\
        $annotation \\
        $args \\
        -i $bam \\
        -o $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandepth: $VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "PANDEPTH module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.8.0'
    """
    touch ${prefix}.chr.stat.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandepth: $VERSION
    END_VERSIONS
    """
}
