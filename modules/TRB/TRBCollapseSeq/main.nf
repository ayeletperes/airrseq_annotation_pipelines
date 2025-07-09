process TRBCollapseSeq {
    tag "TRBCollapseSeq_${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/collapse_files", mode: 'copy'

    input:
        tuple  val(sample_id), path(annotations)
        path v_reference
        val  num_mutation

    output:
        tuple val(sample_id), path("*_collapsed.tsv"), emit: annotations, optional: true
        tuple val(sample_id), path("*_collapsed.fasta"), emit: collapsed_fasta, optional: true

    script:
		
        """
        Rscript "/mnt/bin/rscripts/TRB/TRBCollapseSeq_copy.R" -m ${annotations} -v ${v_reference} -n ${num_mutation}  
        """
}