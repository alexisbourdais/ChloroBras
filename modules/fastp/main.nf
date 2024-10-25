process fastp {

    label 'process_low'

    tag "${sampleId}"

    publishDir "${params.resultsDir}/Trimming_fastp/", mode: 'copy'

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val("${sampleId}"), path("${sampleId}_R{1,2}_trimfastp.fq.gz"), emit: paired_reads
    path("report_fastp.html")

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
    -o ${sampleId}_R1_trimfastp.fq.gz -O ${sampleId}_R2_trimfastp.fq.gz \
    --thread ${task.cpus} \
    --html report_fastp.html
    """
}