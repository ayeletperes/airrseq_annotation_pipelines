process TRBGenotypeInference {
    tag "TRBGenotypeInference_${sample_id}"

    input:
        tuple val(sample_id), path(annotations)
        path v_reference
        path d_reference
        path j_reference
        val min_consensus_count

    output:
        tuple val(sample_id), path("${sample_id}_v_genotype.tsv"), path("${sample_id}_d_genotype.tsv"), path("${sample_id}_j_genotype.tsv"), path("${sample_id}_v_personal.fasta"), path("${sample_id}_d_personal.fasta"), path("${sample_id}_j_personal.fasta"), emit: genotypes

    script:
		
        """
        Rscript "/mnt/bin/rscripts/TRB/TRBGenotypeInference.R" \
        -i ${sample_id} \
        -m ${annotations} \
        -v ${v_reference} \
        -d ${d_reference} \
        -j ${j_reference} \
        -c ${min_consensus_count}
        """
}