process MakeBlastDb {
    tag "MakeBlastDb_${sample_id}_${prefix}"

	input:
        tuple val(sample_id), path(reference_file)
        val  prefix

	output:
		tuple val(sample_id), path("*.db"), emit: blastdb

    script:
        ddb = "${sample_id}_${prefix}.db"
        """
        sed -e '/^>/! s/[.-]//g' ${reference_file} > ${sample_id}_${prefix}_germline.fasta
        touch ${ddb}
        makeblastdb -parse_seqids -dbtype nucl -in ${sample_id}_${prefix}_germline.fasta -out ${ddb}
        """
}