process TRBCollapseSeq {
    tag "TRBCollapseSeq_${sample_id}"

    input:
        tuple  val(sample_id), path(annotations)
        path v_reference

    output:
        tuple val(sample_id), path("*_collapsed.tsv"), emit: annotations, optional: true
        tuple val(sample_id), path("*_collapsed.fasta"), emit: collapsed_fasta, optional: true

    script:
		
        """
        Rscript "/mnt/bin/rscripts/TRB/TRBCollapseSeq.R" -m ${annotations} \
        -v ${v_reference}
        """
}