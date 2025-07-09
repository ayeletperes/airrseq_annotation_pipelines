process DefineClones {

  tag "DefineClones_${sample_id}"

  input:
    tuple val(sample_id), path(reads)


  output:
    tuple val(sample_id), path("*_clone-pass.tsv"), emit: annotations

  script:
    failed      = params.DefineClones.failed == "true" ? "--failed" : ""
    format      = params.DefineClones.format == "airr" ? "--format airr" : "--format changeo"
    seq_field   = params.DefineClones.seq_field ? "--sf ${params.DefineClones.seq_field}" : ""
    v_field     = params.DefineClones.v_field ? "--vf ${params.DefineClones.v_field}" : ""
    d_field     = params.DefineClones.d_field ? "--df ${params.DefineClones.d_field}" : ""
    j_field     = params.DefineClones.j_field ? "--jf ${params.DefineClones.j_field}" : ""
    group_fields= params.DefineClones.group_fields ? "--gf ${params.DefineClones.group_fields}" : ""

    mode   = params.DefineClones.mode   == "gene" ? "" : "--mode ${params.DefineClones.mode}"
    norm   = params.DefineClones.norm   == "len"  ? "" : "--norn ${params.DefineClones.norm}"
    act    = params.DefineClones.act    == "set"  ? "" : "--act ${params.DefineClones.act}"
    model  = params.DefineClones.model  == "ham"  ? "" : "--model ${params.DefineClones.model}"
    sym    = params.DefineClones.sym    == "avg"  ? "" : "--sym ${params.DefineClones.sym}"
    link   = params.DefineClones.link   == "single" ? "" : "--link ${params.DefineClones.link}"
    maxmiss= "--maxmiss ${params.DefineClones.maxmiss}"
    dist = "--dist ${params.DefineClones.dist}"

    """
    DefineClones.py \\
      -d ${reads} \\
      ${failed} \\
      ${format} \\
      ${v_field} \\
      ${d_field} \\
      ${j_field} \\
      ${seq_field} \\
      ${group_fields} \\
      ${mode} \\
      ${act} \\
      ${model} \\
      ${dist} \\
      ${norm} \\
      ${sym} \\
      ${link} \\
      ${maxmiss} \\
      --log DF_${sample_id}.log
    """
}
