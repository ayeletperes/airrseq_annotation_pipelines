process HaplotypeInference {
  tag "haplotype_inference_${sample_id}"

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', saveAs: { fn -> fn.contains('_haplotype.tsv') ? "haplotype/${fn}" : null }
  publishDir "${params.outdir}/${sample_id}", mode: 'copy', saveAs: { fn -> fn.contains('_binomDel.tsv') ? "deletions/${fn}" : null }

  input:
    tuple val(sample_id), path(reads), path(v_germline), path(d_germline)

  output:
    tuple val(sample_id), path("*_haplotype.tsv"), emit: haplotype, optional: true
    tuple val(sample_id), path("*_binomDel.tsv"), emit: deletions, optional: true

  script:

    """
    Rscript /mnt/bin/rscripts/haplotype/rabhit_haplotype_inference.R \\
      --input ${reads} \\
      --v_germline ${v_germline} \\
      --d_germline ${d_germline} \\
      --prefix ${sample_id}
    """
}  
