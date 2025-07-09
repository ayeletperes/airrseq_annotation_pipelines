process MetadataReport {

  publishDir "${params.outdir}", mode: 'copy', saveAs: { fn -> fn.endsWith('.json') ? "meta_data/${fn}" : null }

  input:
    val sample_name
    val aligner               // e.g. "igblast" or "alignair"
    val genotyper             // e.g. "tigger" or "piglet"
    val reference_version     // e.g. "OGRDB revision 8"
    val single_assignment     // "true" or "false"

  output:
    path "annotation_metadata.json", emit: metadata_json

  script:
    """
    Rscript ${workflow.projectDir}/bin/rscripts/new_meta_data.R \\
      --sample_name ${sample_name} \\
      --aligner ${aligner} \\
      --genotyper ${genotyper} \\
      --reference_version '${reference_version}' \\
      --single_assignment ${single_assignment}
    """
}