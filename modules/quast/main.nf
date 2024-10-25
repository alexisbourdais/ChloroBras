process quast {

    label 'process_low'

    publishDir "${params.resultsDir}/", mode: 'copy'

    input:
    val(asm)
    path(allAsm)

    output:
    tuple path("Quast_${asm}/transposed_report.tsv"), path("Quast_${asm}/report.html"), path("Quast_${asm}/report.html"), path("Quast_${asm}/circos/circos.png"), path("Quast_${asm}/genome_stats/"), path("Quast_${asm}/predicted_genes/*.gff")

    script:
    """
    quast -o Quast_${asm} -t ${task.cpus} -r ${params.refFasta} --circos --glimmer --features ${params.refGff} ${allAsm}
    """
}