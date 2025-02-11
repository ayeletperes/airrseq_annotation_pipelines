

params.reads = ""
params.reference_dir = ""

if (!params.reads || !params.reference_dir) {
    throw new IllegalArgumentException("Both 'reads' and 'reference_dir' parameters must be supplied.")
}

params.v_ref = "${params.reference_dir}V_gapped.fasta"
params.d_ref = "${params.reference_dir}D.fasta"
params.j_ref = "${params.reference_dir}J.fasta"
params.aux = "${params.reference_dir}human_gl.aux"
params.ndm = "${params.reference_dir}human.ndm"
sample_name = params.sample_name ?: reads.map { it.name.replaceFirst(/\.[^.]+$/, '') }

include { ReferenceDbWorkflow as InitialReferenceDbWorkflow; ReferenceDbWorkflow as PersonalReferenceDbWorkflow; } from '../../modules/ReferenceDB/workflow/main.nf'
include { IgBlastn as InitialIgBlastn; IgBlastn as PersonalIgBlastn; } from '../../modules/Igblastn/main.nf'
include { MakeDB as InitialMakeDb; MakeDB as PersonalMakeDb; } from '../../modules/MakeDB/main.nf'
include { TRBCollapseSeq } from '../../modules/TRB/TRBCollapseSeq/main.nf'
include { TRBGenotypeInference } from '../../modules/TRB/TRBGenotypeInference/main.nf'
include { TRBGenotypeReport} from '../../modules/TRB/TRBGenotypeReport/main.nf'
include { OGRDBStatsReport } from '../../modules/ogrdb/main.nf'

workflow {
    Channel.fromPath(params.v_ref, type: 'file').set { v_reference }
    Channel.fromPath(params.d_ref, type: 'file').set { d_reference }
    Channel.fromPath(params.j_ref, type: 'file').set { j_reference }
    v_reference.combine(d_reference).combine(j_reference).set { reference }
    InitialReferenceDbWorkflow(v_reference, d_reference, j_reference).set { reference_db }

    
    Channel.fromPath(params.reads, type: 'file').set { reads }
    Channel.fromPath(params.aux, type: 'file').set { aux }
    Channel.fromPath(params.ndm, type: 'file').set { ndm }
    
    // initial Reference run
    InitialIgBlastn(reads, reference_db.db_v_path, reference_db.db_d_path, reference_db.db_j_path, aux, ndm, '_first').set { igblast_out }    
    InitialMakeDb(reads, igblast_out.output, reference, '_first').set { makedb_out }
    TRBCollapseSeq(makedb_out.annotations, v_reference).set { collapse_out }

    // genotype inference
    TRBGenotypeInference(collapse_out.annotations, v_reference, d_reference, j_reference, params.genotype.min_consensus_count).set { genotype_out }
    TRBGenotypeReport(makedb_out.annotations, genotype_out.v_genotypes, genotype_out.d_genotypes, genotype_out.j_genotypes)
    
    genotype_out.v_reference.combine(genotype_out.d_reference).combine(genotype_out.j_reference).set { personal_reference }
    // personal Reference run
    PersonalReferenceDbWorkflow(genotype_out.v_reference, genotype_out.d_reference, genotype_out.j_reference).set { personal_reference_db }
    PersonalIgBlastn(collapse_out.collapsed_reads, personal_reference_db.db_v_path, personal_reference_db.db_d_path, personal_reference_db.db_j_path, aux, ndm, '_personal').set { personal_igblast_out }
    PersonalMakeDb(collapse_out.collapsed_reads, personal_igblast_out.output, personal_reference, '_personal').set { personal_makedb_out }

    // ogrdb reports
    OGRDBStatsReport(personal_makedb_out.annotations, personal_reference)
}