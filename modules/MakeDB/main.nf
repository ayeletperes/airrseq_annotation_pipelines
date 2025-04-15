process MakeDB {
    tag "MakeDB${alignment_suffix}_${sample_id}"

    //publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_MD.*$/) "logs/$filename"}
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', saveAs: { fn -> 
        publish && fn.contains(alignment_suffix) ? "annotations/${fn}" : null
    }
    
    input:
        tuple val(sample_id), path(fasta), path(igblast_out), path(reference)
        val alignment_suffix
        val publish

    output:
        tuple val(sample_id), path("*_db-pass.tsv"), emit: annotations
		path "*_db-fail.tsv", emit: failed, optional: true
		path "*_MD.log", emit: log_file	

    script:
        outname = "${sample_id}${alignment_suffix}"
        ref_files = reference instanceof List ? reference.join(' ') : reference

        failed_opt        = params.makedb.failed         == "true" ? "--failed" : ""
        format_opt        = params.makedb.format         == "changeo" ? "--format changeo" : ""
        extended_opt      = params.makedb.extended       == "true" ? "--extended" : ""
        regions_opt       = params.makedb.regions        == "rhesus-igl" ? "--regions rhesus-igl" : ""
        asisid_opt        = params.makedb.asisid         == "true" ? "--asis-id" : ""
        asiscalls_opt     = params.makedb.asiscalls      == "true" ? "--asis-calls" : ""
        inferjunction_opt = params.makedb.inferjunction  == "true" ? "--infer-junction" : ""
        partial_opt       = params.makedb.partial        == "true" ? "--partial" : ""

        """
        MakeDb.py igblast \\
            -s ${fasta} \\
            -i ${igblast_out} \\
            -r ${ref_files} \\
            --log ${outname}_MD.log \\
            --outname ${outname} \\
            ${extended_opt} \\
            ${failed_opt} \\
            ${format_opt} \\
            ${regions_opt} \\
            ${asisid_opt} \\
            ${asiscalls_opt} \\
            ${inferjunction_opt} \\
            ${partial_opt}
        """
}