process IgBlastn {
    tag "IgBlastn${alignment_suffix}_${sample_id}"

    input:  
        tuple val(sample_id), path(fasta), path(db_v_path), path(db_d_path), path(db_j_path)
        path auxiliary_data
        path custom_internal_data
        val alignment_suffix

    output:
        tuple val(sample_id), path(outfile), emit: output

    script:
    	num_threads = params.igblast.num_threads
		outfmt = params.igblast.outfmt
		num_alignments_V = params.igblast.num_alignments_V
		domain_system = params.igblast.domain_system
        ig_seqtype = params.igblast.ig_seqtype
        d_penalty = params.igblast.d_penalty ? "-D_penalty ${params.igblast.d_penalty}" : ""

		outfile = (outfmt=="MakeDb") ? sample_id + alignment_suffix + ".out" : sample_id + alignment_suffix + ".tsv"
		outfmt = (outfmt=="MakeDb") ? "'7 std qseq sseq btop'" : outfmt

        """
        export IGDATA=/usr/local/share/igblast

        igblastn -query ${fasta} \
            -germline_db_V ${db_v_path.toRealPath()} \
            -germline_db_D ${db_d_path.toRealPath()} \
            -germline_db_J ${db_j_path.toRealPath()} \
            -num_alignments_V ${num_alignments_V} \
            -domain_system ${domain_system} \
            -ig_seqtype ${ig_seqtype} \
            -auxiliary_data ${auxiliary_data} \
            -custom_internal_data ${custom_internal_data} \
            -outfmt ${outfmt} \
            -num_threads ${num_threads} \
            -out ${outfile} \
            ${d_penalty} \
        """
}