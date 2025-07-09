process PIgLETGenotypeInference {

  tag "PIgLET_Genotype_${segment}_${sample_id}"

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', saveAs: { fn -> fn.contains("_genotype_report.tsv") ? "genotype_report/${sample_id}_${segment_flag}_genotype_report.tsv" : null}

  input:
    tuple val(sample_id), path(reads), path(germline)
    path threshold_file
    val segment     
    val find_unmutated
    val single_assignments
    val z_threshold      

  output:
    tuple val(sample_id), path("*_genotype_report.tsv"), emit: genotype
    tuple val(sample_id), path("*_personal_reference.fasta"), emit: reference

  script:
    segment_flag = "${segment.toLowerCase()}_call"

    """
        Rscript /mnt/bin/rscripts/genotype/piglet_genotype_inference.R \\
        --input ${reads} \\
        --germline ${germline} \\
        --call ${segment_flag} \\
        --find_unmutated ${find_unmutated} \\
        --single_assignments ${single_assignments}
    """
}
