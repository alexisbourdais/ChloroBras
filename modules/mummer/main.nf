process mummer {

    label 'process_low'

    publishDir "${results}/Mummer", mode: 'copy'

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    path("${sampleId}_${assembler}.png")

    script:
    """
    nucmer -p ${sampleId} ${params.refFasta} ${assembly}

    mummerplot \
    -x ${params.mummerAxe} \
    -p ${sampleId}_${assembler} \
    -t ${params.mummerFormatOut} \
    "${sampleId}.delta"
    """
}