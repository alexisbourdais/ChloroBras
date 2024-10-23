process getorganelle_index {

    input:
    val index

    output:
    val("")

    script:
    """
    get_organelle_config.py -a ${index}
    """
}

process getorganelle {

    label 'process_high'

    tag "${sampleId}"

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    publishDir "${results}/Assembly/", mode: 'copy', pattern: "*getGood.fasta"

    input:
    val index
    tuple val(sampleId), path(reads)

    output :
    path("${sampleId}_getGood.fasta"), emit: mafft
    tuple val(sampleId), val("get"), path("${sampleId}_get.fasta"), emit: mummer

    script:
    """
    get_organelle_from_reads.py \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${sampleId} \
    -k ${params.getKmer} \
    -F ${params.getIndex} \
    -t ${task.cpus}

    asm1="\$(ls ${sampleId}/embplant_pt.K*.complete.graph1.1.path_sequence.fasta)"
    asm2="\$(ls ${sampleId}/embplant_pt.K*.complete.graph1.2.path_sequence.fasta)"

    rename_fasta_header.py -i "\${asm1}" -n "${sampleId}_get_1" -o "${sampleId}_get_1.fasta"
    rename_fasta_header.py -i "\${asm2}" -n "${sampleId}_get_2" -o "${sampleId}_get_2.fasta"
    cat "${sampleId}_get_1.fasta" "${sampleId}_get_2.fasta" > "${sampleId}_get.fasta"
    selection_assembly.sh "${sampleId}_get_1.fasta" "${sampleId}_get_2.fasta" "${sampleId}_getGood.fasta"
    """
}