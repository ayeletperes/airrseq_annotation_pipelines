process GenerateGappedAlignment {
    tag "GenerateGappedAlignment${alignment_suffix}_${sample_name}"

    input:
    path reads
    path v_germline
    path d_germline
    path j_germline
    val alignment_suffix

    output:
    path output, emit: annotations
    path "*_gapped.log", emit: log_file, optional: true

    script:
    sample_name = params.sample_name ?: reads.getBaseName()
    out_prefix = "${sample_name}${alignment_suffix}"
    output = "${out_prefix}_gapped-passed.tsv"
    is_light = params.chain == "IGH" ? "" : "--is_light"

    """
    python "/mnt/bin/python/igblast_utils/generate_gapped_alignment.py" \\
        --input $reads \\
        --output $output \\
        --v_germline $v_germline \\
        --d_germline $d_germline \\
        --j_germline $j_germline \\
        --log_file ${out_prefix}_gapped.log \\
        ${is_light}
    """
}
