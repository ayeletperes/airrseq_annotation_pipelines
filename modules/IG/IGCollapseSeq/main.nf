process IGCollapseSeq {

    tag "collapse_seq_${sample_id}"

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("*_collapsed-passed.tsv"), emit: annotations
        tuple val(sample_id), path("*_collapsed-seq.fasta"), emit: collapsed_reads
        path "*_collapsed-failed.tsv", emit: failed, optional: true
        path "*_collapsed-duplicate.tsv", emit: duplicates, optional: true
        path "*_collapsed-undetermined.tsv", emit: undetermined, optional: true
        path log_file, emit: log_file, optional: true

    script:
        n_max           = params.Collapse_Seq.n_max ?: 10
        uniq_fields     = params.Collapse_Seq.uniq_fields ?: ['v_call', 'j_call', 'c_call']
        copy_fields     = params.Collapse_Seq.copy_fields ?: ['duplicate_count', 'consensus_count']
        copy_actions    = params.Collapse_Seq.copy_actions ?: ['sum', 'sum']
        maxf            = params.Collapse_Seq.maxf==""? "":"--maxf ${params.Collapse_Seq.maxf} "
        out_prefix      = "${sample_id}"
        log_file        = "${out_prefix}_collapsed.log"
        

        """
        python /mnt/bin/python/IG/CollapseSeq.py \\
            --input ${reads} \\
            --output ${out_prefix} \\
            --max-missing ${n_max} \\
            --uf ${uniq_fields.join(' ')} \\
            --cf ${copy_fields.join(' ')} \\
            --act ${copy_actions.join(' ')} \\
            --log ${log_file} \\
            ${maxf}
        """
}
