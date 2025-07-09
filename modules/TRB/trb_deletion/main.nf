process trbDeletion {
    tag "trbDeletion_${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/deletion", mode: 'copy'

    input:
        tuple val(sample_id), path(annotations)
        path v_reference
        path v_genotype
        path gene_usages_file
        val min_consensus_count

    output:
        tuple val(sample_id), path("${sample_id}_v_genotype_deletion.tsv"), emit: deletion

    script:
		
        """
        Rscript "/mnt/bin/rscripts/TRB/trb_deletion.R" \
        -i ${sample_id} \
        -m ${annotations} \
        -a ${v_reference} \
        -b ${v_genotype} \
        -u ${gene_usages_file} \
        -c ${min_consensus_count}
        """
}
