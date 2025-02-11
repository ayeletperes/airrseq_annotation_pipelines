process IgBlastn {
    input:
        path reads
        path db_v_path
        path db_d_path
        path db_j_path
        path auxiliary_data
        path custom_internal_data
        val alignment_suffix

    output:
        path '*.out', emit: output

    script:
        name = (params.sample_name=="") ? reads.getBaseName() : params.sample_name
    	num_threads = params.igblast.num_threads
		outfmt = params.igblast.outfmt
		num_alignments_V = params.igblast.num_alignments_V
		domain_system = params.igblast.domain_system
        ig_seqtype = params.igblast.ig_seqtype
        d_penalty = params.igblast.d_penalty

		outfile = (outfmt=="MakeDb") ? name + alignment_suffix + ".out" : name + "_" + alignment_suffix + ".tsv"
		outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : outfmt

        """
        export IGDATA=/usr/local/share/igblast

        igblastn -query ${reads} \
            -germline_db_V ${db_v_path.toRealPath()} \
            -germline_db_D ${db_d_path.toRealPath()} \
            -germline_db_J ${db_j_path.toRealPath()} \
            -num_alignments_V ${num_alignments_V} \
            -D_penalty ${d_penalty} \
            -domain_system ${domain_system} \
            -ig_seqtype ${ig_seqtype} \
            -auxiliary_data ${auxiliary_data} \
            -custom_internal_data ${custom_internal_data} \
            -outfmt ${outfmt} \
            -num_threads ${num_threads} \
            -out ${outfile}
        """
}