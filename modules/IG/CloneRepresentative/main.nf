process CloneRepresentative {

  tag "CloneRepresentative_${sample_id}"

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path(outfile), emit: annotations
    path logfile, emit: log, optional: true

  script:
    outfile = "${sample_id}_clone_rep-passed.tsv"
    logfile = "CC_${sample_id}.log"

    """
    Rscript /mnt/bin/rscripts/clones/select_clonal_representatives.R \\
      --input ${reads} \\
      --output ${outfile} \\
      --log ${logfile}
    """
}