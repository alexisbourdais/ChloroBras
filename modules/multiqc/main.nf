process multiqc {

    label 'process_low'

    publishDir "${params.resultsDir}/", mode: 'copy'

    input:
    path(allFastq)

    output :
    path("multiqc_report.html")

    script:
    """
    multiqc ${allFastq}
    """
}