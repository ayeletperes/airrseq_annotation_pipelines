process TRBGenotypeReport {
    publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_genotype.*$/) "genotype/$filename"}

    input:
        path reads
        path v_genotype
        path d_genotype
        path j_genotype

    output:
        path "*_genotype.tsv", emit: genotypes optional true

    script:
        name = params.sample_name
        """
        Rscript "/mnt/bin/rscripts/TRB/TRBGenotypeReport.R" -m ${reads} \
        -v ${v_genotype} \
        -d ${d_genotype} \
        -j ${j_genotype} \
        -o ${name}
        """
}