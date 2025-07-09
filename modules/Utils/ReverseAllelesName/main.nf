process ReverseAllelesName {
    
    tag "reverse_rename_${sample_id}"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', saveAs: { fn -> 
        publish && fn.contains(suffix) ? "rearrangements/${sample_id}${suffix}.tsv" : null
    }

    input:
        tuple val(sample_id), path(annotation_file), path(changes_csv), path(germline_file)
        val suffix
        val region
        val publish

    output:
        tuple val(sample_id), path("${sample_id}${suffix}_${region}_restored_annotation.tsv"), emit: annotations
        tuple val(sample_id), path("${sample_id}${suffix}_${region}_restored_germline.fasta"), emit: germline

    script:       
        """
        cp ${annotation_file} ${sample_id}${suffix}_${region}_restored_annotation.tsv
        cp ${germline_file} ${sample_id}${suffix}_${region}_restored_germline.fasta

        bash /mnt/bin/bash/igblast_utils/reverse_allele_name.sh \\
            ${changes_csv} \\
            ${sample_id}${suffix}_${region}_restored_annotation.tsv \\
            ${sample_id}${suffix}_${region}_restored_germline.fasta
        """
} 
