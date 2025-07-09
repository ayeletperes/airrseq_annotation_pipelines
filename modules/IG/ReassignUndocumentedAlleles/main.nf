process ReassignUndocumentedAlleles {

  tag "Reassign_${call_column}_${sample_id}"

  input:
    tuple val(sample_id), path(reads), path(germline), path(novel)
    val call_column
    val likelihood_column

  output:
    tuple val(sample_id), path("*_reassigned-pass.tsv"), emit: annotations

  script:
    
    """
        python /mnt/bin/python/IG/ReassignAlleles.py \\
        --input ${reads} \\
        --out_name ${sample_id} \\
        --germline_file ${germline} \\
        --novel_file ${novel} \\
        --call_column ${call_column} \\
        --likelihood_column ${likelihood_column}
    """
}
