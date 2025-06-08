process MERQURYFK_PLOIDYPLOT {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/merquryfk:1.1.2--h71df26d_0':
        'biocontainers/merquryfk:1.1.2--h71df26d_0' }"

    input:
    tuple val(meta), path(fastk_hist), path(fastk_ktab)

    output:
    tuple val(meta), path("*.fi.png"), emit: filled_ploidy_plot_png , optional: true
    tuple val(meta), path("*.fi.pdf"), emit: filled_ploidy_plot_pdf , optional: true
    tuple val(meta), path("*.ln.png"), emit: line_ploidy_plot_png   , optional: true
    tuple val(meta), path("*.ln.pdf"), emit: line_ploidy_plot_pdf   , optional: true
    tuple val(meta), path("*.st.png"), emit: stacked_ploidy_plot_png, optional: true
    tuple val(meta), path("*.st.pdf"), emit: stacked_ploidy_plot_pdf, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_PLOIDYPLOT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1.0--h71df26d_1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '1.1.2--h71df26d_0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    PloidyPlot \\
        $args \\
        -T$task.cpus \\
        -o$prefix \\
        ${fastk_ktab.find{ it.toString().endsWith(".ktab") }}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
        r: \$( R --version | sed '1!d; s/.*version //; s/ .*//' )
    END_VERSIONS
    """
}
