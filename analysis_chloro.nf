#!/local/env/env nextflow

/*
===============================================================
 Plastome Analysis Pipeline. Started April 2022.
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

    GetOrganelle
    --getorganelle_index            Index of GetOrganelle, default: "embplant_pt"
    --getorganelle_kmer             Size of kmers, default: "21,45,65,85,105"

    Rename_headers
    --rename_dir_script             Absolute path to script rename

    Nucmer
    --nucmer_ref                    Absolute path to Fasta reference for alignment, default: brassica oleracea

    Mummer
    --mummer_axe                    Size of X-axis (fonction of genome's size), default (plastome): "'[0:154000]'"
    --mummer_format_output          Format of the plot, default: "png"

    Select_assembly
    --select_assembly_script        Absolute path to script select_assembly

    Mafft
    --mafft_method                  Alignment methods, default: "auto"

    Rename_headers_phylo
    --rename_headers_phylo_script   Absolute path to script rename

    Raxml
    --raxml_model                   Model uses by RAxML, default: "GTRGAMMAI"

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
    .set { read_pairs_ch }


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

    errorStrategy 'ignore'

    input:
    val(x) from record_getorganelle_index_ch
    tuple val(pairID), file(reads) from read_pairs_ch

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

process rename_headers {

    tag "${pairID}"

    publishDir "${results}/${pairID}/", mode: 'copy'

    input:
    tuple val(pairID), file(assembly_getorganelle_1), file(assembly_getorganelle_2) from record_getorganelle_ch
    
    output:
    tuple val(pairID), file("${pairID}_get_1"), file("${pairID}_get_2") into record_rename_ch1, record_rename_ch2
    
    script:
    """
    python ${params.rename_script} -i ${assembly_getorganelle_1} -n "${pairID}_get_1" -o "${pairID}_get_1"
    python ${params.rename_script} -i ${assembly_getorganelle_2} -n "${pairID}_get_2" -o "${pairID}_get_2"
    """
}

process concaten {

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_getorganelle_1), file(assembly_getorganelle_2) from record_rename_ch1

    output:
    tuple val(pairID), file("${pairID}.fasta") into record_concat_ch

    script:
    """
    cat ${assembly_getorganelle_1} ${assembly_getorganelle_2} > ${pairID}.fasta
    """
}

process nucmer {

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_cat) from record_concat_ch

    output:
    tuple val(pairID), file("${pairID}.delta") into record_nucmer_ch

    script:
    """
    nucmer -p ${pairID} ${nucmer_ref} ${assembly_cat}
    """
}

process mummer {

    tag "${pairID}"

    publishDir "${results}/${pairID}/", mode: 'move'

    input:
    tuple val(pairID), file(align_from_nucmer) from record_nucmer_ch

    output:
    tuple val(pairID), file("${pairID}.png") into record_mummer_ch

    script:
    """
    mummerplot \
    -x ${params.mummer_axe} \
    -p ${pairID} \
    -t ${params.mummer_format_output} \
    ${align_from_nucmer}
    """
}

process select_assembly {

    tag "${pairID}"

    input:
    tuple val(pairID), file(assembly_getorganelle_1), file(assembly_getorganelle_2) from record_rename_ch2

    output:
    file("${pairID}_graph_right") into record_select_assembly_ch

    script:
    """
    bash ${params.select_assembly_script} ${assembly_getorganelle_1} ${assembly_getorganelle_2} ${pairID}_graph_right
    """
}

process mafft {

    publishDir "${results}/MAFFT", mode: 'copy'

    conda 'mafft'

    tag "Alignment"

    input:
    file(multi_fasta) from record_select_assembly_ch.collectFile(name: 'multi_fasta', newLine: true)

    output:
    file("multi_fasta_align.fasta") into record_mafft_ch

    script:
    """
    mafft --${params.mafft_method} ${multi_fasta} > multi_fasta_align.fasta
    """
}

process rename_headers_phylo {

    tag "Modification of headers"

    input:
    file(multi_fasta_align) from record_mafft_ch

    publishDir "${results}/MAFFT/", mode: 'copy'

    output:
    file("multi_fasta_align_new.fasta") into record_rename_headers_phylo_ch1, record_rename_headers_phylo_ch2

    script:
    """
    bash ${params.rename_headers_phylo_script} ${multi_fasta_align} multi_fasta_align_new.fasta
    """ 
}

process raxml {

    conda 'raxml'

    tag "Creation of the phylogenetic tree"

    publishDir "${results}/RAxML/", mode: 'move'

    input:
    file(multi_fasta_align) from record_rename_headers_phylo_ch1

    output:
    file("RAxML_bootstrap.model_${params.raxml_model}") into record_raxml_ch

    script:
    """
    raxmlHPC \
    -s ${multi_fasta_align} \
    -n model_${params.raxml_model} \
    -d \
    -m ${params.raxml_model} \
    --auto-prot=bic \
    --bootstop-perms=1000 \
    -x 12345 \
    -p 12345 \
    -# autoMRE \
    -f a
    """
}

workflow.onComplete{println("Workflow execution completed sucessfully!")}
