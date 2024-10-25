process mafft {

    label 'process_high'

    publishDir "${params.resultsDir}/Mafft/", mode: 'copy'

    input:
    path(multiFasta)

    output:
    path("align.fasta")

    script:
    """
    mafft --${params.mafftMethod} ${multiFasta} > align.fasta
    sed -i -e 's/_get_1//g' -e 's/_get_2//g' align.fasta
    """
}