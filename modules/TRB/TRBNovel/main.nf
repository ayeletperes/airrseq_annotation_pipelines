process TRBNovel {
    tag "TRBNovel_${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/novel_report", mode: 'copy'

    input:
        tuple  val(sample_id), path(annotations)
        path v_reference

    output:
        tuple val(sample_id), path("igdiscover_novel_selected_igdiscover.tsv"), emit: novel_report, optional: true
        tuple val(sample_id), path("*_novel_V_ref.fasta"), emit: novel_fasta, optional: true


    script:
		
        """
        Rscript "/mnt/bin/rscripts/TRB/TRBnovel.R" -m ${annotations} -v ${v_reference}   
        """
}