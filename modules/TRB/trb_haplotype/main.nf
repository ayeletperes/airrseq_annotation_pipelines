process trbHaplotype {
    tag "trbHaplotype_${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/haplotype", mode: 'copy'

    input:
        tuple val(sample_id), path(annotations)
        path v_reference
        path d_reference
        path genos

    output:
        tuple val(sample_id), path("*_haplotype.tsv"), emit: haplotype, optional: true

    script:
		
        """
        Rscript "/mnt/bin/rscripts/TRB/trb_haplotype.R" \
        --prefix ${sample_id} \
        --input ${annotations} \
        --v_germline ${v_reference} \
        --d_germline ${d_reference} \
        --genos ${genos}
        """
}
