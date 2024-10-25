process mfannot {

    label 'process_medium'

    cache false

    publishDir "${params.resultsDir}/Annotation/", mode: 'move'

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    path("${sampleId}_${assembler}_mfannot.txt")

    script:
    """
    mfannot -g 11 -o ${sampleId}_${assembler}_mfannot.txt ${assembly}
    """
}