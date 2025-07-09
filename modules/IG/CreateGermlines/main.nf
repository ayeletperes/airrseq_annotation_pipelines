process CreateGermlines {

  tag "Create_Germline_${sample_id}"

  input:
    tuple val(sample_id), path(reads), path(reference)
    val cloned

  output:
    tuple val(sample_id), path("*_germ-pass.tsv"), emit: annotations

  script:
    ref_files = reference instanceof List ? reference.join(' ') : reference

    failed      = params.Create_Germlines.failed == "true" ? "--failed" : ""
    format      = params.Create_Germlines.format == "airr" ? "" : "--format changeo"
    g           = "-g ${params.Create_Germlines.g}"
    cloned      = cloned ? "--cloned" : ""
    v_field     = params.Create_Germlines.v_field ? "--vf ${params.Create_Germlines.v_field}" : ""
    d_field     = params.Create_Germlines.d_field ? "--df ${params.Create_Germlines.d_field}" : ""
    j_field     = params.Create_Germlines.j_field ? "--jf ${params.Create_Germlines.j_field}" : ""
    seq_field   = params.Create_Germlines.seq_field ? "--sf ${params.Create_Germlines.seq_field}" : ""
    clone_field = params.Create_Germlines.clone_field ? " ${params.Create_Germlines.clone_field}" : ""

    """
    CreateGermlines.py \\
      -d ${reads} \\
      -r ${ref_files} \\
      ${failed} \\
      ${format} \\
      ${g} \\
      ${cloned} \\
      ${v_field} \\
      ${d_field} \\
      ${j_field} \\
      ${seq_field} \\
      ${clone_field} \\
      --log CG_${sample_id}.log
    """
}
