# Human TRB Annotation Pipeline

## Overview

The **Human TRB Annotation Pipeline** is a Nextflow-based workflow designed for **T cell receptor beta (TRB) annotation and genotype inference**. It processes **FASTA files**, infers genotypes, and OGRDB stats reports.

---

## Required Docker / Singularity Images

Before running the pipeline, pull the necessary containers:

```bash
docker pull peresay/base:latest
docker pull williamlees/ogrdbstats:latest
```

Note: it is not necesseray to pre-pull the singularity images, the pipeline will automatically pull them

```bash
singularity pull docker://peresay/base:latest
singularity pull docker://williamlees/ogrdbstats:latest
```

---

## Execution Profiles

The pipeline supports multiple execution profiles:

### `docker` Profile

- Uses **Docker** to manage dependencies and runs on **local** executer
- Requires to manually pull the docker containers

### `singularity` Profile

- Uses **Singularity** instead of Docker and runs on **local** executer
- Automatically pulls containers if not available

### `singularityHPC` Profile

- Runs on **SLURM** clusters with Singularity
- Automatically pulls containers if not available

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
    profile=singularity # Specify the execution profile
    nextflow -log $log_file run ../airrseq_annotation_pipelines/pipelines/Human_TRB_Annotation/main.nf \
        -w $work_dir \
        --reference_dir $reference_dir \
        --sample_name $sample_name \
        --reads $reads \
        -profile $profile
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

## Configuration Parameters

The pipeline supports configurable parameters, defined in `nextflow.config`. Below are the key parameters:

### General Parameters

- **`outdir`**: Output directory (default: `results`)
- **`ncors`**: Number of CPU threads to use (default: `4`)

### IgBlastn Parameters

- **`num_threads`**: Number of CPU threads to use (default: `params.ncors`)
- **`outfmt`**: Output format (default: `MakeDb`)
- **`num_alignments_V`**: Maximum V calls (default: `10`)
- **`domain_system`**: Domain annotation system (default: `imgt`)
- **`ig_seqtype`**: Sequence type (default: `TCR`)
- **`d_penalty`**: D alignment penalty (default: `-2`)

### MakeDb Parameters

- **`failed`**: Output failed sequences (`false` by default)
- **`format`**: Output format (`airr` by default)
- **`regions`**: Annotation regions (`default` by default)
- **`extended`**: Use extended annotation (`true` by default)
- **`asisid`**: Use original sequence ID (`false` by default)
- **`asiscalls`**: Use original calls (`false` by default)
- **`inferjunction`**: Infer junctions (`false` by default)
- **`partial`**: Include partial alignment results (`false` by default)

### Genotype Inference Parameters

- **`min_consensus_count`**: Minimum consensus count for sequence inclusion (`1` by default)

### OGRDB Stats Parameters

- **`chain`**: TCR chain type (`TRBV` by default)

---
