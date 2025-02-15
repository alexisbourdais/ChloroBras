process mummer {

    label 'process_low'

    publishDir "${params.resultsDir}/Mummer", mode: 'copy'

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    path("${sampleId}_${assembler}.png")

    script:
    """
    nucmer -p ${sampleId} ${params.refFasta} ${assembly}

    mummerplot \
    -x '[0:154000]' \
    -p ${sampleId}_${assembler} \
    -t "png" \
    "${sampleId}.delta"
    """
}
