# Human TRB Annotation Pipeline

## Overview

The **Human TRB Annotation Pipeline** is a Nextflow-based workflow designed for **T cell receptor beta (TRB) annotation and genotype inference**. It processes **FASTA files**, infers genotypes, and OGRDB stats reports.

---

## Required Docker Images

Before running the pipeline, pull the necessary containers:

```bash
docker pull peresay/base:latest
docker pull williamlees/ogrdbstats:latest
```

---

## Running the Pipeline

### **1. Organize Input Data**

Ensure your **FASTA input files** are in the correct directory. Example:

```plaintext
project_root/
├── airrseq_annotation_pipelines/
│   ├── pipelines/
│   │   ├── TRB_Annotation_Pipeline/
│   │   │   ├── main.nf
│   │   │   ├── nextflow.config
│   ├── modules/
│   ├── bin/
├── input/
│   ├── C4.fasta
```

### **2. Run the Pipeline with a Single File**

1. Create a directory for the work and results

    ```bash
    sample_dir="C4"
    mkdir -p $sample_dir
    ```

2. Move to the working directory to ensure cache location

    ```bash
    cd $sample_dir
    ```

3. Run analysis

    ```bash
    reads="../C4.fasta" # AIRR-seq file
    sample_name="C4" # Sample name
    log_file=logs/nextflow.log # log file location. Optional
    reference_dir="../TRB/rev1/" # Path to reference directory
    work_dir=work # Processing directory name
    nextflow -log $log_file run ../airrseq_annotation_pipelines/pipelines/Human_TRB_Annotation/main.nf \
        -w $work_dir \
        --reference_dir $reference_dir \
        --sample_name $sample_name \
        --reads $reads \
        -with-docker # this ensure the pipeline will run with the dockers
    ```

---

## Expected Output

Results are saved in the `results/` directory and include:

| File Path                                      | Description                                      |
|------------------------------------------------|--------------------------------------------------|
| `annotations/*_personal.tsv`                   | Annotated reads with personal reference set      |
| `genotype/_genotype.tsv`                       | Genotype inference results                       |
| `ogrdbstats/*_ogrdb_report.tsv`                | OGRDB statistics reports                         |
| `ogrdbstats/*_ogrdb_plots.tsv`                 | OGRDB statistics plots                           |

---

## Debugging & Logs

To check logs:

```bash
nextflow log test/nextflow.log
```

For debugging input files:

```nextflow
println "Params.reads: ${params.reads}"
println "Processed Reads: ${reads}"
```

---
