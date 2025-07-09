process InferUndocumentedAlleles {

    tag "infer_novel_${sample_id}"

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', saveAs: { fn -> fn.contains("_novel-passed.tsv") ? "novel_report/${sample_id}_novel-passed.tsv" : null}

    input:
        tuple val(sample_id), path(reads), path(v_germline_file)

    output:
        tuple val(sample_id), path("${sample_id}_novel-passed.tsv"), emit: novel_table, optional: true
        tuple val(sample_id), path("V_novel_germline.fasta"), emit: novel_fasta, optional: true
        tuple val(sample_id), env('NOVEL_FOUND'), emit: novel_found, optional: true
        
    script:
        script_path = params.Undocumented_Alleles.chain == 'IGH' ? '/mnt/bin/rscripts/novel_alleles/IGHNovelAllelesInference.R' : '/mnt/bin/rscripts/novel_alleles/IGLNovelAllelesInference.R'
        num_threads   = params.Undocumented_Alleles.num_threads
        germline_min  = params.Undocumented_Alleles.germline_min
        min_seqs      = params.Undocumented_Alleles.min_seqs
        y_intercept   = params.Undocumented_Alleles.y_intercept
        alpha         = params.Undocumented_Alleles.alpha
        j_max         = params.Undocumented_Alleles.j_max
        min_frac      = params.Undocumented_Alleles.min_frac
        auto_mutrange = params.Undocumented_Alleles.auto_mutrange
        mut_range     = params.Undocumented_Alleles.mut_range
        pos_range     = params.Undocumented_Alleles.pos_range
        cmd           = params.Undocumented_Alleles.undocumented_alleles_use_custom_tigger ? 'source("/mnt/bin/rscripts/novel_alleles/functions_tigger.R")' : 'library(tigger)'

        if(params.skip_novel_inference){
            """
            export NOVEL_FOUND="no_novel"
            """
        }else{
            """
            Rscript ${script_path} \
                ${reads} \
                ${v_germline_file} \
                ${sample_id}_novel-passed.tsv \
                V_novel_germline \
                '${cmd}' \
                ${num_threads} \
                ${germline_min} \
                ${min_seqs} \
                ${y_intercept} \
                ${alpha} \
                ${j_max} \
                ${min_frac} \
                ${auto_mutrange} \
                ${mut_range} \
                ${pos_range}
        
            if [ -s "${sample_id}_novel-passed.tsv" ]; then
                export NOVEL_FOUND="novel"
            else
                export NOVEL_FOUND="no_novel"
            fi
            """
        }
        
}
