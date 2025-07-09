include { IGCollapseSeq }                                    from '../IGCollapseSeq/main.nf'
include { IgBlastAlignmentWorkflow as InitialIgBlastAlignmentWorkflow } from '../IgBlastAlignmentWorkflow/main.nf'
include { IgBlastAlignmentWorkflow as NovelIgBlastAlignmentWorkflow   } from '../IgBlastAlignmentWorkflow/main.nf'
include { InferUndocumentedAlleles }                          from '../InferUndocumentedAlleles/main.nf'
include { IGFilterAnnotations }                               from '../IGFilterAnnotations/main.nf'

workflow IgBlastAlignmentAndNovelWorkflow {

    /*───────────── inputs ─────────────*/
    take:
        alignment_input   // (sid, reads, v_ref, d_ref, j_ref)               – FASTA per sample
        aux               // aux file (single path)
        ndm               // ndm file (single path)
        publish_initial   // publishDir bool
        publish_novel     // publishDir bool

    /*───────────── main DAG ───────────*/
    main:

        /* 1 ─ initial IgBlast alignment */
        InitialIgBlastAlignmentWorkflow(
            alignment_input,          // (sid, reads)
            aux,
            ndm,
            '_first',
            publish_initial
        ).set { init_align_raw }           // (sid, ann, v_ref, d_ref, j_ref)

        /* split annotations vs. references, keep key */
        init_align_raw.map { sid, ann, _v, _d, _j -> tuple(sid, ann) }
                      .set { ann_init }     // (sid, ann)

        init_align_raw.map { sid, _ann, v, d, j -> tuple(sid, v, d, j) }
                      .set { ref_init }     // (sid, v, d, j)

        /* 2 ─ collapse & filter reads    */
        IGCollapseSeq(ann_init)                        // two keyed sub‑channels (sid, ann‑collapsed), (sid, fasta‑collapsed)
            .set { collapsed }

        IGFilterAnnotations(collapsed.annotations)     // (sid, ann‑collapsed)
            .set { filt }                              // filt.annotations  (sid, ann‑filt)
        filt.filtered_reads.set { filtered_reads }     // (sid, reads‑filt)
        
        /* 3 ─ undocumented allele inference */
        // build (sid, reads‑filt, v_ref) input
        filt.annotations
            .join( ref_init.map{ sid, v, _d, _j -> tuple(sid, v) } )
            .map { sid, ann, v -> tuple(sid, ann, v) }
            .set { undoc_input }

        InferUndocumentedAlleles(undoc_input)
            .set { novel_out }                         // ─► novel_out.novel_fasta   (sid, fasta?) opt
                                                       //    novel_out.novel_found   (sid, flag)

        /* 4 ─ always run “novel” alignment            */
        
        filtered_reads
            .join(novel_out.novel_fasta)
            .join( ref_init.map{ sid, _v, d, j -> tuple(sid, d, j) } )
            .map { sid, reads, novel_fasta, d, j -> tuple(sid, reads, novel_fasta, d, j) }
            .set { novel_input }

        NovelIgBlastAlignmentWorkflow(
            novel_input, // (sid, reads, novel_found, v_ref, d_ref, j_ref)
            aux,
            ndm,
            '_novel',
            publish_novel
        ).set { novel_align }                         // (sid, ann, v_ref, d_ref, j_ref)

        /* 5 ─ decide per‑sample which path to keep    */
        novel_out.novel_found                         // (sid, "novel" | "no_novel")
            .branch { _sid, flag ->
                USE_NOVEL   : flag == 'novel'
                USE_FALLBACK: flag == 'no_novel'
            }
            .set { novel_switch }

        /* 5a ─ choose annotation table */
        novel_switch.USE_NOVEL
            .join( novel_align.map{ sid, ann, _v, _d, _j -> tuple(sid, ann) } )
            .map { sid, _flag, ann -> tuple(sid, ann) }
            .set { ann_novel }

        novel_switch.USE_FALLBACK
            .join( filt.annotations )
            .map { sid, _flag, ann -> tuple(sid, ann) }
            .set { ann_fallback }

        ann_novel.mix(ann_fallback)
                 .set { clone_annotations }           // (sid, ann‑final)

        /* 5b ─ choose germline reference set */
        novel_switch.USE_NOVEL
            .join( novel_align.map{ sid, _ann, v, d, j -> tuple(sid, v, d, j) } )
            .map { sid, _flag, v, d, j -> tuple(sid, v, d, j) }
            .set { ref_novel }

        novel_switch.USE_FALLBACK
            .join( ref_init )
            .map { sid, _flag, v, d, j -> tuple(sid, v, d, j) }
            .set { ref_fallback }

        ref_novel.mix(ref_fallback)
                 .set { clone_reference }             // (sid, v, d, j)

        /* ─ optional debug prints ─ */
        clone_annotations.view { it -> "Annotation ➜ ${it[0]} : ${it[1]}" }
        clone_reference  .view { it -> "References  ➜ ${it[0]} : ${it[1..3]}" }

    /*───────────── outputs ───────────*/
    emit:
        filtered_reads    // (sid, reads‑filt)
        clone_annotations // (sid, ann‑final)
        clone_reference   // (sid, v, d, j)
}