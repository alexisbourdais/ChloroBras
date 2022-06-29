#!/local/env/env nextflow

/*
===============================================================
 Plastome Analysis Pipeline. Started April 2022.
 Testing with Nextflow v21.10.6
 #### Homepage / Documentation
 https://github.com/Airgetlam35/ChloroBras/
 #### Authors
 Alexis Bourdais
---------------------------------------------------------------
*/

def helpMessage() {
    log.info """

    Defines the pipeline inputs parameters
    Each of the following parameters can be specified as command line options or in the config file
    Beware of quotes when necessary


    Reads directory
    --readsFiles                    Path to input data, default: "./samples/*_R{1,2}.fastq.gz"

    Results directory
    --resultsDir                    Path to results directory, default: "./results/"

    Sqtk
    --seqtk_nb_read                 Number of reads to keep, default: 2000000

    GetOrganelle
    --getorganelle_index            Index of GetOrganelle, default: "embplant_pt"
    --getorganelle_kmer             Size of kmers, default: "21,45,65,85,105"

    FastPlast
    --fastplast_index               Index of Fast-Plast, default: "Brassicales"

    OrgAsm
    --orgasm_probes                 Index of ORGanelle ASeMbler, default: "protChloroArabidopsis"

    Rename_headers
    --rename_dir_script             Absolute path to script rename

    Nucmer
    --nucmer_ref                    Absolute path to Fasta reference for alignment, default: brassica oleracea

    Mummer
    --mummer_axe                    Size of X-axis (fonction of genome's size), default (plastome): "'[0:154000]'"
    --mummer_format_output          Format of the plot, default: "png"

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
        helpMessage()
        exit 0
}


///////////////////////////////////////////////////////////
//////////////////////    Variables    ////////////////////
///////////////////////////////////////////////////////////


results = file(params.resultsDir)
nucmer_ref = file(params.nucmer_ref)

Channel
    .fromFilePairs(params.readsFiles)
    .ifEmpty {error "Cannot find any reads matching: ${params.readsFiles}"}
    .into { read_pairs_ch1; read_pairs_ch2 }


///////////////////////////////////////////////////////////
//////////////////////    Processes    ////////////////////
///////////////////////////////////////////////////////////

process getorganelle_index {
    
    tag "${params.getorganelle_index}"

    label 'getorganelle'

    output:
    val("") into record_getorganelle_index_ch

    script:
    """
    get_organelle_config.py --add ${params.getorganelle_index}
    """
}

process getorganelle {

    label 'getorganelle'

    tag "${pairID}"

    time '2h'

    errorStrategy 'ignore'

    input:
    val(x) from record_getorganelle_index_ch
    tuple val(pairID), file(reads) from read_pairs_ch1

    output:
    tuple val(pairID), file("${pairID}/embplant_pt.K*.complete.graph1.1.path_sequence.fasta"), file("${pairID}/embplant_pt.K*.complete.graph1.2.path_sequence.fasta") optional true into record_getorganelle_ch

    script:
    """
    get_organelle_from_reads.py \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${pairID} \
    -k ${params.getorganelle_kmer} \
    -F ${params.getorganelle_index} \
    -t ${task.cpus}
    """
}

process rename_headers_get {

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_getorganelle_1), file(assembly_getorganelle_2) from record_getorganelle_ch
    
    output:
    tuple val(pairID), file("${pairID}_get_1"), file("${pairID}_get_2") into record_rename_get_ch
    
    script:
    """
    python ${params.rename_script} -i ${assembly_getorganelle_1} -n "${pairID}_get_1" -o "${pairID}_get_1"
    python ${params.rename_script} -i ${assembly_getorganelle_2} -n "${pairID}_get_2" -o "${pairID}_get_2"
    """
}

process concaten {

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_getorganelle_1), file(assembly_getorganelle_2) from record_rename_get_ch

    output:
    tuple val(pairID), file("${pairID}.fasta") into record_concat_get_ch

    script:
    """
    cat ${assembly_getorganelle_1} ${assembly_getorganelle_2} > ${pairID}.fasta
    """
}

process nucmer_get {

    label 'nucmer'

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_cat) from record_concat_get_ch

    output:
    tuple val(pairID), file("${pairID}.delta") into record_nucmer_get_ch

    script:
    """
    nucmer -p ${pairID} ${nucmer_ref} ${assembly_cat}
    """
}

process mummer_get {

    label 'mummer'

    tag "${pairID}"

    publishDir "${results}/${pairID}/", mode: 'move'

    input:
    tuple val(pairID), file(align_from_nucmer) from record_nucmer_get_ch

    output:
    tuple val(pairID), file("${pairID}_get.png") into record_mummer_get_ch

    script:
    """
    mummerplot \
    -x ${params.mummer_axe} \
    -p ${pairID}_get \
    -t ${params.mummer_format_output} \
    ${align_from_nucmer}
    """
}

process seqtk {

    tag "${pairID}"
    
    input:
    tuple val(pairID), file(reads) from read_pairs_ch2

    output:
    tuple val(pairID), file("${pairID}_sub1.fq.gz"), file("${pairID}_sub2.fq.gz") into record_seqtk_ch1, record_seqtk_ch2

    script:
    """
    seqtk sample -s666 ${reads[0]} ${params.seqtk_nb_read} |gzip -c - > ${pairID}_sub1.fq.gz
    seqtk sample -s666 ${reads[1]} ${params.seqtk_nb_read} |gzip -c - > ${pairID}_sub2.fq.gz
    """
}

process orgasm_index {

    label 'orgasm'

    tag "${pairID}"

    input:
    tuple val(pairID), file(sample_R1), file(sample_R2) from record_seqtk_ch1

    output:
    tuple val(pairID), file("${pairID}.odx") into record_orgasm_index_ch

    script:
    """
    oa index --estimate-length=0.9 ${pairID} ${sample_R1} ${sample_R2}
    """
}

process orgasm_buildgraph {

    label 'orgasm'

    tag "${pairID}"

    time '2h'

    errorStrategy 'ignore'

    input:
    tuple val(pairID), file(index) from record_orgasm_index_ch

    output:
    tuple val(pairID), file(index), file("${pairID}.oas") optional true into record_orgasm_build_ch

    script:
    """
    oa buildgraph --probes ${params.orgasm_probes} ${pairID} ${pairID}
    """
}

process orgasm_unfold {

    label 'orgasm'

    tag "${pairID}"

    errorStrategy 'ignore'

    input:
    tuple val(pairID), file(index), file(graph) from record_orgasm_build_ch

    output:
    tuple val(pairID), file("${pairID}_org.fasta") into record_orgasm_unfold_ch

    script:
    """
    oa unfold ${pairID} ${pairID} > ${pairID}_org.fasta
    """
}

process rename_headers_org {

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_orgasm) from record_orgasm_unfold_ch
    
    output:
    tuple val(pairID), file("${pairID}_org") into record_rename_org_ch
    
    script:
    """
    python ${params.rename_script} -i ${assembly_orgasm} -n "${pairID}_org" -o "${pairID}_org"
    """
}

process nucmer_org {

    label 'nucmer'

    tag "${pairID}"

    errorStrategy 'ignore'

    input:
    tuple val(pairID), file(assembly_org) from record_rename_org_ch

    output:
    tuple val(pairID), file("${pairID}.delta") into record_nucmer_org_ch

    script:
    """
    nucmer -p ${pairID} ${nucmer_ref} ${assembly_org}
    """
}

process mummer_org {

    label 'mummer'

    tag "${pairID}"

    publishDir "${results}/${pairID}/", mode: 'move'

    input:
    tuple val(pairID), file(align_from_nucmer) from record_nucmer_org_ch

    output:
    tuple val(pairID), file("${pairID}_org.png") into record_mummer_org_ch

    script:
    """
    mummerplot \
    -x ${params.mummer_axe} \
    -p ${pairID}_org \
    -t ${params.mummer_format_output} \
    ${align_from_nucmer}
    """
}

process fastplast {

    tag "${pairID}"

    time '2h'

    errorStrategy 'ignore' 

    input:
    tuple val(pairID), file(sample_R1), file(sample_R2) from record_seqtk_ch2

    output:
    tuple val(pairID), file("${pairID}/Final_Assembly/${pairID}_FULLCP.fsa") optional true into record_fastplast_ch

    script:
    """
    fast-plast.pl \
    -1 ${sample_R1} \
    -2 ${sample_R2} \
    --name ${pairID} \
    --bowtie_index ${params.fastplast_index} \
    """
}

process rename_headers_fast {

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_fastplast) from record_fastplast_ch
    
    output:
    tuple val(pairID), file("${pairID}_fast") into record_rename_fast_ch
    
    script:
    """
    python ${params.rename_script} -i ${assembly_fastplast} -n "${pairID}_fast" -o "${pairID}_fast"
    """
}

process nucmer_fast {

    label 'nucmer'

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_cat) from record_rename_fast_ch

    output:
    tuple val(pairID), file("${pairID}.delta") into record_nucmer_fast_ch

    script:
    """
    nucmer -p ${pairID} ${nucmer_ref} ${assembly_cat}
    """
}

process mummer_fast {

    label 'mummer'

    tag "${pairID}"

    publishDir "${results}/${pairID}/", mode: 'move'

    input:
    tuple val(pairID), file(align_from_nucmer) from record_nucmer_fast_ch

    output:
    tuple val(pairID), file("${pairID}_fast.png") into record_mummer_fast_ch

    script:
    """
    mummerplot \
    -x ${params.mummer_axe} \
    -p ${pairID}_fast \
    -t ${params.mummer_format_output} \
    ${align_from_nucmer}
    """
}

workflow.onComplete{println("Workflow execution completed sucessfully!")}
