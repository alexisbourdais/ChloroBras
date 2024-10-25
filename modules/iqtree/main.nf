process iqtree {

    label 'process_high'

    errorStrategy 'ignore'

    publishDir "${params.resultsDir}/Iqtree/", mode: 'copy'

    input:
    path(multi_fasta_align)

    output:
    path("iqtree_${params.iqtreeModel}.*")

    script:
    if (params.iqtreeOption=="") {
    """
    iqtree \
    -s ${multi_fasta_align} \
    -pre "iqtree_${params.iqtreeModel}" \
    -m ${params.iqtreeModel} \
    -T ${task.cpus}
    """
    }

    else {
    """
    iqtree \
    -s ${multi_fasta_align} \
    -pre "iqtree_${params.iqtreeModel}" \
    -m ${params.iqtreeModel} \
    -T ${task.cpus} \
    ${params.iqtreeOption}
    """
    }
}