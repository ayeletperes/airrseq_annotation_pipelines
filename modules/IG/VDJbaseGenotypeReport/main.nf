process VDJbaseGenotypeReport {
  tag "VDJbaseGenotypeReport_${sample_id}"
  
  publishDir "${params.outdir}/${sample_id}", mode: 'copy', saveAs: { fn -> fn.contains("_Final_genotype.tsv") ? "genotype_report/${sample_id}_Final_genotype.tsv" : null}

  input:
    tuple val(sample_id), path(initial_reads), path(final_reads), path(v_genotype), path(d_genotype), path(j_genotype)

  output:
    tuple val(sample_id), path("*_Final_genotype.tsv"), emit: genotype_report

  when:
    params.generate_vdjbase_genotype_report

  script:

    """
    Rscript /mnt/bin/rscripts/vdjbase_reports/vdjbase_genotype_report.R \\
      --initial_reads ${initial_reads} \\
      --personal_reads ${final_reads} \\
      --v_genotype ${v_genotype} \\
      --d_genotype ${d_genotype} \\
      --j_genotype ${j_genotype} \\
      --outname ${sample_id}
    """
}  
