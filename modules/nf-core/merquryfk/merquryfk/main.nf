process MERQURYFK_MERQURYFK {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/merquryfk:1.1.2--h71df26d_0':
        'biocontainers/merquryfk:1.1.2--h71df26d_0' }"

    input:
    tuple val(meta), path(fastk_hist), path(fastk_ktab), path(fastk_data), path(assembly), path(haplotigs)
    tuple val(met),  path(mat_ktab), path(mat_data)                                                         //optional
    tuple val(met),  path(pat_ktab), path(pat_data)                                                         //optional

    output:
    tuple val(meta), path("${prefix}.completeness.stats")         , emit: stats
    tuple val(meta), path("${prefix}.*_only.bed")                 , emit: bed
    tuple val(meta), path("${prefix}.*.qv")                       , emit: assembly_qv
    tuple val(meta), path("${prefix}.*.spectra-cn.fl.{png,pdf}")  , emit: spectra_cn_fl,  optional: true
    tuple val(meta), path("${prefix}.*.spectra-cn.ln.{png,pdf}")  , emit: spectra_cn_ln,  optional: true
    tuple val(meta), path("${prefix}.*.spectra-cn.st.{png,pdf}")  , emit: spectra_cn_st,  optional: true
    tuple val(meta), path("${prefix}.qv")                         , emit: qv
    tuple val(meta), path("${prefix}.spectra-asm.fl.{png,pdf}")   , emit: spectra_asm_fl, optional: true
    tuple val(meta), path("${prefix}.spectra-asm.ln.{png,pdf}")   , emit: spectra_asm_ln, optional: true
    tuple val(meta), path("${prefix}.spectra-asm.st.{png,pdf}")   , emit: spectra_asm_st, optional: true
    tuple val(meta), path("${prefix}.phased_block.bed")           , emit: phased_block_bed,   optional: true
    tuple val(meta), path("${prefix}.phased_block.stats")         , emit: phased_block_stats, optional: true
    tuple val(meta), path("${prefix}.continuity.N.{pdf,png}")     , emit: continuity_N,       optional: true
    tuple val(meta), path("${prefix}.block.N.{pdf,png}")          , emit: block_N,            optional: true
    tuple val(meta), path("${prefix}.block.blob.{pdf,png}")       , emit: block_blob,         optional: true
    tuple val(meta), path("${prefix}.hapmers.blob.{pdf,png}")     , emit: hapmers_blob,       optional: true
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_MERQURYFK module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1.0--h71df26d_1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.2--h71df26d_0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    mkdir -p tmp
    
    MerquryFK \\
        $args \\
        -T$task.cpus \\
        -Ptmp \\
        $fastk_ktab \\
        $mat_ktab \\
        $pat_ktab \\
        $assembly \\
        $haplotigs \\
        $prefix

    rm -r tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1.0--h71df26d_1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.1--h71df26d_1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.completeness.stats
    touch ${prefix}.qv
    touch ${prefix}._.qv
    touch ${prefix}._only.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
