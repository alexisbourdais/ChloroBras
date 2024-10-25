process raxmlng {

    label 'process_high'

    errorStrategy 'ignore'

    publishDir "${params.resultsDir}/RAxML-NG/", mode: 'copy'

    input:
    path(multi_fasta_align)

    output:
    path("raxmlng_${params.raxmlngModel}.raxml.*")

    script:

    if (params.raxmlngOption=="") {
    """
    raxml-ng \
    --all \
    --bs-trees ${params.raxmlngBootstrap} \
    --msa ${multi_fasta_align} \
    --model ${params.raxmlngModel} \
    --prefix "raxmlng_${params.raxmlngModel}" \
    --threads ${task.cpus}
    """
    }

    else {
    """
    raxml-ng \
    --all \
    --bs-trees ${params.raxmlngBootstrap} \
    --msa ${multi_fasta_align} \
    --model ${params.raxmlngModel} \
    --prefix "raxmlng_${params.raxmlngModel}" \
    --threads ${task.cpus} \
    ${params.raxmlngOption==""}
    """
    }
}