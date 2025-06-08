process SMUDGEPLOT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${projectDir}/assets/containers/fastk_merquryfk_smudgeplot.sif"

    input:
    tuple val(meta), path(ktab), path(data)
    val(cutoff)

    output:
    tuple val(meta), path("*.smu")            , emit: smu
    tuple val(meta), path("*.pdf")            , emit: pdf
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p tmp
    
    smudgeplot.py \\
        hetmers \\
        -t $task.cpus \\
        -tmp tmp \\
        -L $cutoff \\
        -o $prefix \\
        $args \\
        $ktab

    smudgeplot.py \\
        all \\
        -o $prefix \\
        $args2 \\
        ${prefix}_text.smu

    rm -r tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smudgeplot: \$(smudgeplot.py --version 2>&1 | sed 's/^.*smudgeplot //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_text.smu
    touch ${prefix}_centralities.pdf
    touch ${prefix}_smudgeplot.pdf
    touch ${prefix}_smudgeplot_log10.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        smudgeplot: \$(smudgeplot.py --version 2>&1 | sed 's/^.*smudgeplot //')
    END_VERSIONS
    """
}
