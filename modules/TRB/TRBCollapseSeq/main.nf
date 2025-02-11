process TRBCollapseSeq {
    input:
        path reads
        path v_reference

    output:
        path "*_collapsed.tsv", emit: annotations optional true
        path "*_collapsed.fasta", emit: collapsed_reads optional true

    script:
		
        """
        Rscript "/mnt/bin/rscripts/TRB/TRBCollapseSeq.R" -m ${reads} \
        -v ${v_reference}
        """
}