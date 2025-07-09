process IGFilterAnnotations {

    tag "filter_annotations_${sample_id}"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("*_filtered-passed.tsv"), emit: annotations
        tuple val(sample_id), path("*_filtered-seq.fasta"), emit: filtered_reads, optional: true
        path "*_filtered-failed.tsv", emit: failed, optional: true
        path log_file, emit: log_file, optional: true

    script:
        n_max           = params.Filter_Annot.n_max ?: 10
        min_consensus   = params.Filter_Annot.min_consensus ?: 0
        fasta           = params.Filter_Annot.fasta ?: " --fasta"

        out_prefix      = "${sample_id}"
        log_file        = "${out_prefix}_filter.log"

        """
        python /mnt/bin/python/IG/IGFilterAnnotations.py \\
            --input ${reads} \\
            --output ${out_prefix} \\
            --n-limit ${n_max} \\
            --min-consensus ${min_consensus} \\
            ${fasta} \\
            --log ${log_file}
        """
}
