/*───────────────────────────────────────────────────────────────*
 *  AlignAIRAlignmentAndNovelWorkflow (multi‑sample version)     *
 *───────────────────────────────────────────────────────────────*/

include { IGCollapseSeq             } from '../IGCollapseSeq/main.nf'
include { AlignAIR                  } from '../../AlignAIR/main.nf'
include { InferUndocumentedAlleles  } from '../InferUndocumentedAlleles/main.nf'
include { IGFilterAnnotations       } from '../IGFilterAnnotations/main.nf'
include { ReassignUndocumentedAlleles } from '../ReassignUndocumentedAlleles/main.nf'

workflow AlignAIRAlignmentAndNovelWorkflow {

    /*───────────── input ─────────────*/
    take:
        alignment_input        // (sid, reads, v_ref, d_ref, j_ref)

    /*───────────── main DAG ──────────*/
    main:

        /* 1 ─ AlignAIR first pass */
        alignair_in = alignment_input.map { sid, reads, _v, _d, _j ->
                        tuple(sid, reads) }                 // (sid, reads)

        AlignAIR( alignair_in, '_first' )
            .set { align_raw }                              // (sid, ann)

        /* 2 ─ collapse + filter */
        IGCollapseSeq( align_raw )
            .set { collapsed }                              // two keyed sub‑channels (sid, ann‑collapsed), (sid, fasta‑collapsed)

        IGFilterAnnotations( collapsed.annotations )        // (sid, ann‑collapsed)
            .set { filt }                                   // filt.annotations  (sid, ann‑filt)

        filt.filtered_reads.set { filtered_reads }                   // filt.filtered_reads (sid, reads‑filt)

        /* 3 ─ undocumented‑allele inference  */
        v_ref_ch = alignment_input.map { sid, _reads, v, _d, _j ->
                        tuple(sid, v) }                     // (sid, v_ref)

        undoc_input = filt.annotations
                        .join( v_ref_ch )                   // key = sid
                        .map { sid, ann, v -> tuple(sid, ann, v) }

        InferUndocumentedAlleles( undoc_input )
            .set { novel_out }                              // novel_out.novel_table (sid, tsv?)
                                                            // novel_out.novel_fasta (sid, fasta?)
                                                            // novel_out.novel_found (sid, flag)

        /* 4 ─ re‑assign calls if novel alleles present     */
        // Build (sid, reads, germline, novel) tuple
        reassign_input = filt.annotations
                            .join( v_ref_ch )               // add original V germline
                            .join( novel_out.novel_table )  // add novel table
                            .map { sid, ann, v, novel_tbl ->
                                   tuple(sid, ann, v, novel_tbl) }

        ReassignUndocumentedAlleles(
            reassign_input,
            'v_call',
            'v_likelihoods'
        ).set { reassign_out }                              // (sid, ann‑reassigned)

        /* 5 ─ decide per‑sample whether to use novel path  */
        novel_out.novel_found                               // (sid, 'novel' | 'no_novel')
            .branch { _sid, flag ->
                USE_NOVEL   : flag == 'novel'
                USE_FALLBACK: flag == 'no_novel'
            }
            .set { novel_switch }

        /* 5a ─ pick annotation table */
        ann_novel = novel_switch.USE_NOVEL
                        .join( reassign_out )               // (sid, flag, ann‑reassigned)
                        .map  { sid, _flag, ann -> tuple(sid, ann) }

        ann_fallback = novel_switch.USE_FALLBACK
                        .join( filt.annotations )
                        .map  { sid, _flag, ann -> tuple(sid, ann) }

        ann_novel.mix(ann_fallback)
                 .set { clone_annotations }                 // (sid, ann‑final)

        /* 5b ─ pick germline reference trio (V,D,J) */
        dj_ref_ch = alignment_input.map { sid, _reads, _v, d, j ->
                        tuple(sid, d, j) }                  // (sid, d_ref, j_ref)

        ref_novel = novel_switch.USE_NOVEL
                        .join( novel_out.novel_fasta )      // add novel V
                        .join( dj_ref_ch )                  // add D & J
                        .map { sid, _flag, v_novel, d, j ->
                               tuple(sid, v_novel, d, j) }

        ref_fallback = novel_switch.USE_FALLBACK
                        .join( alignment_input.map{ sid, _r, v, d, j -> tuple(sid, v, d, j) } )
                        .map { sid, _flag, v, d, j -> tuple(sid, v, d, j) }

        ref_novel.mix(ref_fallback)
                 .set { clone_reference }                   // (sid, v, d, j)

        /* optional debug */
        clone_annotations.view{ it -> "ANN  ➜ ${it[0]} : ${it[1]}" }
        clone_reference  .view{ it -> "REF  ➜ ${it[0]} : ${it[1..3]}" }

    /*───────────── outputs ───────────*/
    emit:
        filtered_reads    // (sid, reads‑filt)
        clone_annotations // (sid, ann‑final)
        clone_reference   // (sid, v, d, j)
}
