process orgasm {

    label 'process_high'

    publishDir "${results}/Assembly/", mode: 'copy'

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    tag "${sampleId}"

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), val("org"), path("${sampleId}_org.fasta"), emit: mummer
    path("${sampleId}_org.fasta"), emit: mafft

    script:
    """
    oa index --estimate-length=0.9 ${sampleId} ${reads[0]} ${reads[1]}
    oa buildgraph --probes ${params.orgasmProbes} ${sampleId} ${sampleId}
    oa unfold ${sampleId} ${sampleId} > ${sampleId}_orgasm.fasta
    convert_multiline_oneline.sh "${sampleId}_orgasm.fasta" "${sampleId}_oneLine_org.fasta"
    rename_fasta_header.py -i "${sampleId}_oneLine_org.fasta" -n "${sampleId}" -o "${sampleId}_org.fasta"
    """
}