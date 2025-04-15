process OGRDBStatsReport {
    tag "OGRDBStatsReport_${sample_id}"

    publishDir "${params.outdir}/${sample_id}/ogrdbstats", mode: 'copy'

    input:
        tuple val(sample_id), path(annotations), path(reference)
        val suffix // output file suffix

    output:
        path "*csv", emit: stats
        path "*pdf", emit: report

    script:
        chain = params.ogrdbstats.chain

        reference = reference instanceof List ? reference.join(' ') : reference
        output_prefix = "${sample_id}${suffix}"

        """
        cat ${reference} > ${output_prefix}_reference_set.fasta
        germline_file_path=\$(realpath ${output_prefix}_reference_set.fasta)

        novel=""

        if grep -q "_[A-Z][0-9]" \${germline_file_path}; then
            awk '/^>/{f=0} \$0 ~ /_[A-Z][0-9]/ {f=1} f' \${germline_file_path} > ${output_prefix}_novel_sequences.fasta
            novel=\$(realpath ${output_prefix}_novel_sequences.fasta)
            diff \$germline_file_path \$novel | grep '^<' | sed 's/^< //' > ${output_prefix}_ogrdbstats_germline.fasta
            germline_file_path=\$(realpath ${output_prefix}_ogrdbstats_germline.fasta)
            novel="--inf_file \$novel"
        fi

        ogrdbstatsReads=${annotations}
        IFS='\t' read -a var < ${annotations}

        if [[ ! "\${var[*]}" =~ "v_call_genotyped" ]]; then
            awk -F'\t' '{col=\$5;gsub("call", "call_genotyped", col); print \$0 "\\t" col}' ${annotations} > ${output_prefix}_genotyped.tsv
            ogrdbstatsReads=${output_prefix}_genotyped.tsv
        fi

        ogrdbstatsReads_path=\$(realpath \$ogrdbstatsReads)

        run_ogrdbstats \\
            \$germline_file_path \\
            "Homosapiens" \\
            \$ogrdbstatsReads_path \\
            ${chain} \\
            \$novel
        """
}