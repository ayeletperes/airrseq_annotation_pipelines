process MakeDB {
    //publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_MD.*$/) "logs/$filename"}
    publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}_personal.*$/) "annotations/$filename"}
    
    input:
        path reads
        path igblast_out
        path reference
        val alignment_suffix

    output:
        path "*_db-pass.tsv", emit: annotations
		path "*_db-fail.tsv", emit: failed optional true
		path "${name}_MD*", emit: log_file	

    script:
        // params
		name = (params.sample_name=="") ? reads.getBaseName() : params.sample_name
		failed = params.makedb.failed
		format = params.makedb.format
		regions = params.makedb.regions
		extended = params.makedb.extended
		asisid = params.makedb.asisid
		asiscalls = params.makedb.asiscalls
		inferjunction = params.makedb.inferjunction
		partial = params.makedb.partial

		failed = (failed=="true") ? "--failed" : ""
		format = (format=="changeo") ? "--format changeo" : ""
		extended = (extended=="true") ? "--extended" : ""
		regions = (regions=="rhesus-igl") ? "--regions rhesus-igl" : ""
		asisid = (asisid=="true") ? "--asis-id" : ""
		asiscalls = (asiscalls=="true") ? "--asis-calls" : ""
		inferjunction = (inferjunction=="true") ? "--infer-junction" : ""
		partial = (partial=="true") ? "--partial" : ""

        // inputs
		outname = name + alignment_suffix
        reference = reference instanceof List ? reference.join(' ') : reference 

        """
        MakeDb.py igblast \
                -s ${reads} \
                -i ${igblast_out} \
                -r ${reference} \
                --log ${name}_MD.log \
                --outname ${outname}\
                ${extended} \
                ${failed} \
                ${format} \
                ${regions} \
                ${asisid} \
                ${asiscalls} \
                ${inferjunction} \
                ${partial}
        """
}