process OGRDBStatsReport {
    publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename.matches(".*_ogrdb_(report|plots).*")) "ogrdbstats/${filename}"}

    input:
        path reads
        path reference

    output:
        path "*csv", emit: stats
        path "*pdf", emit: report

    script:
        name = params.sample_name
        chain = params.ogrdbstats.chain

        reference = reference instanceof List ? reference.join(' ') : reference
        """
        
        cat ${reference} > reference_set.fasta;
        
        germline_file_path=\$(realpath reference_set.fasta)

        novel=""

        if grep -q "_[A-Z][0-9]" \${germline_file_path}; then
            awk '/^>/{f=0} \$0 ~ /_[A-Z][0-9]/ {f=1} f' \${germline_file_path} > novel_sequences.fasta
            novel=\$(realpath novel_sequences.fasta)
            diff \$germline_file_path \$novel | grep '^<' | sed 's/^< //' > ogrdbstats_germline.fasta
            germline_file_path=\$(realpath ogrdbstats_germline.fasta)
            novel="--inf_file \$novel"
        fi

        IFS='\t' read -a var < ${reads}

        ogrdbstatsReads=${reads}

        if [[ ! "\${var[*]}" =~ "v_call_genotyped" ]]; then
            awk -F'\t' '{col=\$5;gsub("call", "call_genotyped", col); print \$0 "\t" col}' ${reads} > ${name}_genotyped.tsv
            ogrdbstatsReads=${name}_genotyped.tsv
        fi

        ogrdbstatsReads_path=\$(realpath \$ogrdbstatsReads)


        run_ogrdbstats \
            \$germline_file_path \
            "Homosapiens" \
            \$ogrdbstatsReads_path \
            ${chain} \
            \$novel 
        """
}