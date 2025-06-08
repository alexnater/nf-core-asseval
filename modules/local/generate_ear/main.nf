process GENERATE_EAR {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta) , path(template)
    tuple val(meta2), path(genomescope), path(smudgeplot)
    tuple val(meta3), path(hifi_depth), path(ul_depth), path(hic_depth)
    tuple val(meta4), path(stats), path(qv)
    tuple val(meta5), path(busco)

    output:
    tuple val(meta), path("${prefix}_ear.yaml"), emit: yaml
    tuple val(meta), path("${prefix}_ear.pdf") , emit: pdf
    path "versions.yml"                        , emit: versions

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    def smudgeplot_arg = smudgeplot ? "--smudgeplot $smudgeplot" : ''

    """
    generate_yaml.py \\
        --genomescope $genomescope \\
        $smudgeplot_arg \\
        --merqury \$PWD \\
        --hifi $hifi_depth \\
        --ul $ul_depth \\
        --hic $hic_depth \\
        $args \\
        -o ${prefix}_ear.yaml \\
        $template
        
    make_EAR.py ${prefix}_ear.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}