process makeNdm {
    tag "makeNdm_${sample_id}_${reference_file}"

  	input:
     tuple val(sample_id), path(reference_file)
  
  	output:
  		path("*.ndm"), emit: ndmfile

    script:
        """
        make_igblast_ndm ${reference_file} VB human.ndm
        """
}