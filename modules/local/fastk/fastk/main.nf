process FASTK_FASTK {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastk:1.1.0--h71df26d_1' :
        'biocontainers/fastk:1.1.0--h71df26d_1' }"

    input:
    tuple val(meta), path(reads)
    val kvalue

    output:
    tuple val(meta), path("*.hist")                      , emit: hist
    tuple val(meta), path("*.txt")                       , emit: txt
    tuple val(meta), path("*.ktab*")                     , emit: ktab, optional: true
    tuple val(meta), path(".*.ktab*", hidden: true)      , emit: data, optional: true
    tuple val(meta), path("*.{prof,pidx}*", hidden: true), emit: prof, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_FASTK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1.0--h71df26d_1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    mkdir -p tmp

    FastK \\
        -k$kvalue \\
        -T$task.cpus \\
        -M${task.memory.toGiga()} \\
        -N${prefix} \\
        -Ptmp \\
        $args \\
        $reads

    rm -r tmp

    Histex \\
        $args2 \\
        ${prefix}.hist \\
        > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_ktab = args.contains('-t') ? "touch ${prefix}_fk.ktab .${prefix}_fk.ktab.1" : ''
    def touch_prof = args.contains('-p') ? "touch ${prefix}_fk.prof .${prefix}_fk.pidx.1" : ''
    """
    touch ${prefix}.hist
    touch ${prefix}.txt
    $touch_ktab
    $touch_prof

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """
}
