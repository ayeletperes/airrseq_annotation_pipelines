process TRBGenotypeReport {
    tag "TRBGenotypeReport_${sample_id}"

    publishDir "${params.outdir}/${sample_id}/genotype", mode: 'copy'
    
    input:
        tuple val(sample_id), path(annotations), path(v_genotype), path(d_genotype), path(j_genotype)

    output:
        path "*_genotype.tsv", emit: genotypes, optional: true

    script:

        """
        Rscript "/mnt/bin/rscripts/TRB/TRBGenotypeReport.R" -m ${annotations} \
        -v ${v_genotype} \
        -d ${d_genotype} \
        -j ${j_genotype} \
        -o ${sample_id}
        """
}