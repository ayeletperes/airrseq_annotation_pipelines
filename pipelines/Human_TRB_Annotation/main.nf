include { MakeBlastDb as MakeBlastDbV; MakeBlastDb as MakeBlastDbD; MakeBlastDb as MakeBlastDbJ; MakeBlastDb as MakeBlastDbPersonalV; MakeBlastDb as MakeBlastDbPersonalD; MakeBlastDb as MakeBlastDbPersonalJ } from '../../modules/MakeBlastDB/main.nf'
include { ShortenAllelesName as V_ShortenAllelesName; ShortenAllelesName as D_ShortenAllelesName; ShortenAllelesName as J_ShortenAllelesName} from '../../modules/Utils/ShortenAllelesName/main.nf'
include { ReverseAllelesName as V_ReverseAllelesName; ReverseAllelesName as D_ReverseAllelesName; ReverseAllelesName as J_ReverseAllelesName} from '../../modules/Utils/ReverseAllelesName/main.nf'
include { IgBlastn as InitialIgBlastn; IgBlastn as PersonalIgBlastn; } from '../../modules/Igblastn/main.nf'
include { MakeDB as InitialMakeDb; MakeDB as PersonalMakeDb; } from '../../modules/MakeDB/main.nf'
include { TRBCollapseSeq } from '../../modules/TRB/TRBCollapseSeq/main.nf'
include { TRBGenotypeInference } from '../../modules/TRB/TRBGenotypeInference/main.nf'
include { TRBGenotypeReport } from '../../modules/TRB/TRBGenotypeReport/main.nf'
include { OGRDBStatsReport } from '../../modules/ogrdb/main.nf'
include { makeNdm } from '../../modules/makeNdm/main.nf'
include { IgBlastAlignmentWorkflow as InitialAlignment; IgBlastAlignmentWorkflow as PersonalAlignment} from '../../modules/IG/IgBlastAlignmentWorkflow/main.nf'

// Assign reference paths
params.reference = [
    v_ref: "${params.reference_dir}/V_gapped.fasta",
    d_ref: "${params.reference_dir}/D.fasta",
    j_ref: "${params.reference_dir}/J.fasta",
    aux: "${params.reference_dir}/human_gl.aux",
    ndm: "${params.reference_dir}/human.ndm"
]

// Wrap everything inside this:
workflow {
    // Load reference files once
    v_reference = file(params.reference.v_ref)
    d_reference = file(params.reference.d_ref)
    j_reference = file(params.reference.j_ref)
    aux_file = file(params.reference.aux)
    ndm_file = file(params.reference.ndm)
    // Create a channel of samples
    if (params.samples) {
        // Multi-sample mode: read from samples file
        sample_files = Channel
            .fromPath(params.samples, checkIfExists: true)
            .map { tuple(it.baseName, file(it)) }
    } else {
        if (params.input_csv) {
            // Multi-sample mode: read from samples file
            sample_files = Channel
                .fromPath(params.input_csv, checkIfExists: true)
                .splitCsv(header: true, sep: ',')
                .map { row -> tuple(row.sample_id, file(row.file)) }
        } else {
            // Single sample mode
            sample_files = Channel
                .fromPath(params.fasta, checkIfExists: true)
                .map { tuple(it.baseName, file(it)) }
        }
    }

    sample_files
            .map { it -> tuple(it[0], it[1], v_reference, d_reference, j_reference) }
            .set { alignment_input }

    suffix = "_first"

    InitialAlignment(
            alignment_input, // (sid, reads, v_ref, d_ref, j_ref)
            aux_file,
            ndm_file,
            "_first",
            false).set { init_align_raw }

    init_align_raw.map { sid, ann, _v, _d, _j -> tuple(sid, ann) }
                      .set { ann_init }     // (sid, ann)

    init_align_raw.map { sid, _ann, v, d, j -> tuple(sid, v, d, j) }
                      .set { ref_init }     // (sid, v, d, j)

    // Run Collapse
    TRBCollapseSeq(ann_init, v_reference, 3).set { collapse_out_raw }


    // Check for low depth samples if enabled
    if(params.check_low_depth) {
        collapse_out_raw.annotations
            .map { sample_id, annotations_file ->
                def count = annotations_file.countLines() - 1  // Subtract 1 for header
                tuple(sample_id, annotations_file, count)
            }
            .branch {
                pass: it[2] >= params.check_low_depth_threshold
                fail: it[2] < params.check_low_depth_threshold
            }
            .set { depth_check }
        
        // Log warning for discarded samples
        depth_check.fail
            .map { sample_id, _unused, count ->
                log.warn "⚠️ Sample ${sample_id} discarded due to low depth (${count} < ${params.check_low_depth_threshold})"
                sample_id
            }
            .set { _discarded_samples }
        
        // Continue with samples that passed the threshold
        depth_check.pass
            .map { sample_id, annotations_file, _unused -> tuple(sample_id, annotations_file) }
            .set { collapse_out_annotations }
        
        // Combine with collapsed_fasta, keeping only samples that passed
        collapse_out_raw.collapsed_fasta
            .join(collapse_out_annotations, by: 0)
            .map { sample_id, fasta, _annotations -> tuple(sample_id, fasta) }
            .set { collapse_out_collapsed_fasta }
    } else {        // If check_low_depth is disabled, use all samples
        collapse_out_annotations = collapse_out_raw.annotations
        collapse_out_collapsed_fasta = collapse_out_raw.collapsed_fasta
    }


    // Run Genotype
    TRBGenotypeInference(
        collapse_out_annotations,
        v_reference,
        d_reference,
        j_reference,
        params.genotype.min_consensus_count,
        params.genotype.max_snp_position
    ).set { genotype_out }

    // Create Genotype Report
    genotype_out.map { it -> tuple(it[0], it[1], it[2], it[3]) }.set { genotypes }
    ann_init.join(genotypes).set{ genotype_report_input }
    TRBGenotypeReport(genotype_report_input)

    // Create personal reference
    genotype_out.map { it -> tuple(it[0], it[4]) }.set { v_reference_ch }
    genotype_out.map { it -> tuple(it[0], it[5]) }.set { d_reference_ch }
    genotype_out.map { it -> tuple(it[0], it[6]) }.set { j_reference_ch }

    // Now pass the four channels to your workflow:
    personal_db_v_path = MakeBlastDbPersonalV(v_reference_ch, 'V')
    personal_db_d_path = MakeBlastDbPersonalD(d_reference_ch, 'D')
    personal_db_j_path = MakeBlastDbPersonalJ(j_reference_ch, 'J')

    collapse_out_collapsed_fasta
    .join(personal_db_v_path, by: 0)
    .join(personal_db_d_path, by: 0)
    .join(personal_db_j_path, by: 0)
    .map { sample_id, sample_file, v_db, d_db, j_db ->
        tuple(sample_id, sample_file, v_db, d_db, j_db)
    }
    .set { personal_igblast_input }

    // Run IgBlast
    PersonalIgBlastn(
        personal_igblast_input,
        aux_file,
        ndm_file,
        '_personal'
    ).set { personal_igblast_raw }

    // Run MakeDB

    collapse_out_collapsed_fasta.join(personal_igblast_raw)
    .join(v_reference_ch).join(d_reference_ch).join(j_reference_ch)
        .map { sample_id, fasta, annotations, v_file, d_file, j_file ->
            tuple(sample_id, fasta, annotations, [v_file, d_file, j_file])
        }
        .set { personal_makedb_input }

    PersonalMakeDb(personal_makedb_input, '_personal', true).set { personal_makedb_out }

    // Create OGRDB report
    // personal_makedb_out.annotations
    // .join(v_reference_ch).join(d_reference_ch).join(j_reference_ch)
    //     .map { sample_id, annotations, v_file, d_file, j_file ->
    //         tuple(sample_id, annotations, [v_file, d_file, j_file])
    //     }
    //     .set { personal_ogrdb_input }

    personal_makedb_out.annotations
        .map { it -> tuple(it[0], it[1], [v_reference, d_reference, j_reference]) }
        .set { personal_ogrdb_input }


    OGRDBStatsReport(personal_ogrdb_input, "")
}