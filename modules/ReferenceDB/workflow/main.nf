include { MakeBlastDb as MakeBlastDbV; MakeBlastDb as MakeBlastDbD; MakeBlastDb as MakeBlastDbJ } from '../MakeBlastDB/main.nf'

workflow ReferenceDbWorkflow {
    take:
        v_reference
        d_reference
        j_reference

    main:
        db_v_path = MakeBlastDbV(v_reference, true)
        db_d_path = MakeBlastDbD(d_reference, db_v_path.ready)
        db_j_path = MakeBlastDbJ(j_reference, db_d_path.ready)

    emit:
        db_v_path = db_v_path.blastdb
        db_d_path = db_d_path.blastdb
        db_j_path = db_j_path.blastdb
}