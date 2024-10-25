process trimgalore {

    label 'process_low'

    tag "${sampleId}"

    publishDir "${params.resultsDir}/Trimming_trimgalore/", mode: 'copy'

    input:
    tuple val(sampleId), path(reads)

    output :
    tuple val("${sampleId}"), path("${sampleId}_R{1,2}_trimgalore.fq.gz"), emit: paired_reads
    tuple path("${sampleId}_R1_trimgalore_fastqc.zip"), path("${sampleId}_R2_trimgalore_fastqc.zip"), emit: qc
    tuple path("${sampleId}_R1_trimgalore_report.txt"), path("${sampleId}_R2_trimgalore_report.txt")

    script:
    """
    trim_galore \
    --paired ${reads[0]} ${reads[1]} \
    --basename ${sampleId} \
    --gzip \
    -o "Trimming_trimgalore" \
    --cores "${task.cpus}" \
    --fastqc
    #--illumina --stranded_illumina

    mv "Trimming_trimgalore/${sampleId}_val_1.fq.gz" "${sampleId}_R1_trimgalore.fq.gz"
    mv "Trimming_trimgalore/${sampleId}_val_2.fq.gz" "${sampleId}_R2_trimgalore.fq.gz"
    mv "Trimming_trimgalore/${sampleId}_val_1_fastqc.zip" "${sampleId}_R1_trimgalore_fastqc.zip"
    mv "Trimming_trimgalore/${sampleId}_val_2_fastqc.zip" "${sampleId}_R2_trimgalore_fastqc.zip"
    mv "Trimming_trimgalore/${reads[0]}_trimming_report.txt" "${sampleId}_R1_trimgalore_report.txt"
    mv "Trimming_trimgalore/${reads[1]}_trimming_report.txt" "${sampleId}_R2_trimgalore_report.txt"
    """
}