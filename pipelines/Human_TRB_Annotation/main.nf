include { MakeBlastDb as MakeBlastDbV; MakeBlastDb as MakeBlastDbD; MakeBlastDb as MakeBlastDbJ; MakeBlastDb as MakeBlastDbPersonalV; MakeBlastDb as MakeBlastDbPersonalD; MakeBlastDb as MakeBlastDbPersonalJ } from '../../modules/MakeBlastDB/main.nf'
include { IgBlastn as InitialIgBlastn; IgBlastn as PersonalIgBlastn; } from '../../modules/Igblastn/main.nf'
include { MakeDB as InitialMakeDb; MakeDB as PersonalMakeDb; } from '../../modules/MakeDB/main.nf'
include { TRBCollapseSeq } from '../../modules/TRB/TRBCollapseSeq/main.nf'
include { TRBGenotypeInference } from '../../modules/TRB/TRBGenotypeInference/main.nf'
include { TRBGenotypeReport } from '../../modules/TRB/TRBGenotypeReport/main.nf'
include { OGRDBStatsReport } from '../../modules/ogrdb/main.nf'

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

    // create igblast reference
    v_reference_ch = sample_files.map { it -> tuple(it[0], v_reference) }
    d_reference_ch = sample_files.map { it -> tuple(it[0], d_reference) }
    j_reference_ch = sample_files.map { it -> tuple(it[0], j_reference) }

    db_v_path = MakeBlastDbV(v_reference_ch, 'V')
    db_d_path = MakeBlastDbD(d_reference_ch, 'D')
    db_j_path = MakeBlastDbJ(j_reference_ch, 'J')
    
    sample_files
    .join(db_v_path, by: 0)
    .join(db_d_path, by: 0)
    .join(db_j_path, by: 0)
    .map { sample_id, sample_file, v_db, d_db, j_db ->
        tuple(sample_id, sample_file, v_db, d_db, j_db)
    }
    .set { igblast_input }

    // Run IgBlast
    InitialIgBlastn(
        igblast_input,
        aux_file,
        ndm_file,
        '_first'
    ).set { igblast_out_raw }

    // Run MakeDB
    sample_files.join(igblast_out_raw).map { it -> tuple(it[0], it[1], it[2], [v_reference, d_reference, j_reference]) }.set{ makedb_input }

    InitialMakeDb(makedb_input, '_first', false).set { makedb_out_raw }

    // Run Collapse
    TRBCollapseSeq(makedb_out_raw.annotations, v_reference).set { collapse_out_raw }


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
        params.genotype.min_consensus_count
    ).set { genotype_out }

    // Create Genotype Report
    genotype_out.map { it -> tuple(it[0], it[1], it[2], it[3]) }.set { genotypes }
    makedb_out_raw.annotations.join(genotypes).set{ genotype_report_input }
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
    personal_makedb_out.annotations
    .join(v_reference_ch).join(d_reference_ch).join(j_reference_ch)
        .map { sample_id, annotations, v_file, d_file, j_file ->
            tuple(sample_id, annotations, [v_file, d_file, j_file])
        }
        .set { personal_ogrdb_input }

    OGRDBStatsReport(personal_ogrdb_input, "")
}