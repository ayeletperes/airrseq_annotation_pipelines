process MakeBlastDb {
	input:
		path reference_file
		val ready

	output:
		path "*.db", emit: blastdb
		val true, emit: ready

    script:
        gdb = reference_file.getBaseName() + '.fasta'
		ddb = reference_file.getBaseName() + '.db'
        """
            sed -e '/^>/! s/[.-]//g' ${reference_file} > germline.fasta
            touch ${ddb}
            makeblastdb -parse_seqids -dbtype nucl -in germline.fasta -out ${ddb}
        """
}