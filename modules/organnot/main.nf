process organnot {

    label 'process_medium'

    publishDir "${params.resultsDir}/Annotation/", mode: 'copy'

    input:
    tuple val(sampleId), val(assembler), path(assembly)

    output:
    path("${sampleId}_${assembler}_organnot.txt")

    script:
    if (params.workflow == "fromAsm") {
        """
        organnot -c ${assembly}
        mv "${sampleId}.annot.circular" "${sampleId}_${assembler}_organnot.txt"
        """
    }
    else {
        """
        organnot -c ${assembly}
        mv ${sampleId}_${assembler}.annot.circular ${sampleId}_${assembler}_organnot.txt
        """
    }
}