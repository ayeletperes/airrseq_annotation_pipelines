process ShortenAllelesName {

    tag "rename_fasta_${sample_id}"

    input:
        tuple val(sample_id), path(input_fasta)
        val suffix

    output:
        tuple val(sample_id), path("${sample_id}${suffix}_renamed.fasta"), emit: renamed_fasta
        tuple val(sample_id), path("${sample_id}${suffix}_changes.csv"), emit: changes_log, optional: true

    script:
        
        """
        python3 /mnt/bin/python/igblast_utils/shorten_alleles_names.py \
            ${input_fasta} \
            ${sample_id}${suffix}_renamed.fasta \
            ${sample_id}${suffix}_changes.csv
        """
} 
