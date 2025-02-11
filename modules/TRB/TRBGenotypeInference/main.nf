process TRBGenotypeInference {
    input:
        path reads
        path v_reference
        path d_reference
        path j_reference
        val min_consensus_count

    output:
        path "v_genotype.tsv", emit: v_genotypes optional true
        path "d_genotype.tsv", emit: d_genotypes optional true
        path "j_genotype.tsv", emit: j_genotypes optional true
        path "v_personal.fasta", emit: v_reference optional true
        path "d_personal.fasta", emit: d_reference optional true
        path "j_personal.fasta", emit: j_reference optional true

    script:
		
        """
        Rscript "/mnt/bin/rscripts/TRB/TRBGenotypeInference.R" -m ${reads} \
        -v ${v_reference} \
        -d ${d_reference} \
        -j ${j_reference} \
        -c ${min_consensus_count}
        """
}