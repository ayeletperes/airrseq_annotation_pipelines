include { MakeBlastDb as MakeBlastDbV; MakeBlastDb as MakeBlastDbD; MakeBlastDb as MakeBlastDbJ } from '../../MakeBlastDB/main.nf'
include { IgBlastn } from '../../Igblastn/main.nf'
include { MakeDB } from '../../MakeDB/main.nf'
include { ShortenAllelesName as V_ShortenAllelesName; ShortenAllelesName as D_ShortenAllelesName; ShortenAllelesName as J_ShortenAllelesName} from '../../Utils/ShortenAllelesName/main.nf'
include { ReverseAllelesName as V_ReverseAllelesName; ReverseAllelesName as D_ReverseAllelesName; ReverseAllelesName as J_ReverseAllelesName} from '../../Utils/ReverseAllelesName/main.nf'

workflow IgBlastAlignmentWorkflow {

  take:
    data
    aux_file
    ndm_file
    suffix
    publish

  main:

    // unpack data tuple
    sample_files = data.map { it -> tuple(it[0], it[1]) }

    v_reference = data.map { it -> tuple(it[0], it[2]) }
    d_reference = data.map { it -> tuple(it[0], it[3]) }
    j_reference = data.map { it -> tuple(it[0], it[4]) }
      
    V_ShortenAllelesName(v_reference, suffix+"_V").set { shorten_v_reference }
    D_ShortenAllelesName(d_reference, suffix+"_D").set { shorten_d_reference }
    J_ShortenAllelesName(j_reference, suffix+"_J").set { shorten_j_reference }

    db_v_path = MakeBlastDbV(shorten_v_reference.renamed_fasta, 'V')
    db_d_path = MakeBlastDbD(shorten_d_reference.renamed_fasta, 'D')
    db_j_path = MakeBlastDbJ(shorten_j_reference.renamed_fasta, 'J')
    
    sample_files
    .join(db_v_path, by: 0)
    .join(db_d_path, by: 0)
    .join(db_j_path, by: 0)
    .map { sample_id, sample_file, v_db, d_db, j_db ->
        tuple(sample_id, sample_file, v_db, d_db, j_db)
    }
    .set { igblast_input }

    IgBlastn(
      igblast_input,
      aux_file,
      ndm_file,
      suffix
    ).set { igblast_out_raw }

    sample_files.join(igblast_out_raw)
    .join(shorten_v_reference.renamed_fasta)
    .join(shorten_d_reference.renamed_fasta)
    .join(shorten_j_reference.renamed_fasta)
    .map { it -> tuple(it[0], it[1], it[2], [it[3], it[4], it[5]]) }.set{ makedb_input }

    MakeDB(
      makedb_input,
      '_first',
      publish
    ).set { makedb_out_raw }

    makedb_out_raw.annotations
    .join(shorten_v_reference.changes_log)
    .join(shorten_v_reference.renamed_fasta)
    .set{ v_reverseAllelesName_input }

    V_ReverseAllelesName(
      v_reverseAllelesName_input,
      suffix,
      "V",
      false
      ).set { v_restored_outputs }

    v_restored_outputs.annotations
    .join(shorten_d_reference.changes_log)
    .join(shorten_d_reference.renamed_fasta)
    .set{ d_reverseAllelesName_input }

    D_ReverseAllelesName(
      d_reverseAllelesName_input,
      suffix,
      "D",
      false
      ).set { d_restored_outputs }

    d_restored_outputs.annotations
    .join(shorten_j_reference.changes_log)
    .join(shorten_j_reference.renamed_fasta)
    .set{ j_reverseAllelesName_input }

    J_ReverseAllelesName(
      j_reverseAllelesName_input,
      suffix,
      "J",
      publish
      ).set { j_restored_outputs }
    
    j_restored_outputs.annotations
    .join(v_restored_outputs.germline)
    .join(d_restored_outputs.germline)
    .join(j_restored_outputs.germline)
    .set { alignment_output }

  emit:
    alignment_output

} 
