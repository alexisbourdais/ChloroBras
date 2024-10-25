process fastplast {

    label 'process_high'

    publishDir "${params.resultsDir}/Assembly/", mode: 'copy'

    tag "${sampleId}"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), val("fpt"), path("${sampleId}_fpt.fasta"), emit: mummer
    path("${sampleId}_fpt.fasta"), emit: mafft

    script:
    """
    fast-plast.pl \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    --name ${sampleId} \
    --bowtie_index ${params.fastIndex} \
    --threads ${task.cpus}

    mv "${sampleId}/Final_Assembly/${sampleId}_FULLCP.fsa" "${sampleId}_fpt.fasta"
    """
}