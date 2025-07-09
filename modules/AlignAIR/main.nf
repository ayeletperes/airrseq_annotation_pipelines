process AlignAIR {
    tag "AlignAIR${alignment_suffix}_${sample_id}"

    input:
        tuple val(sample_id), path(reads)
        val alignment_suffix

    output:
        tuple val(sample_id), path("*_alignair_results.tsv"), emit: annotations

    script:
        chain_type = params.chain=="IGH" ? "heavy" : "light"
    	model_checkpoint = chain_type=="heavy"? "/app/pretrained_models/IGH_S5F_576" : "/app/pretrained_models/IGL_S5F_576"

		outfile = sample_id + alignment_suffix
		save_path = "./${outfile}_"

        """
        python /app/app.py run --sequences ${reads} --model-checkpoint ${model_checkpoint} --chain-type ${chain_type} --airr-format --save-path ${save_path}
        """
}