# Human TRB Annotation Pipeline

## Overview

The **Human TRB Annotation Pipeline** is a Nextflow-based workflow designed for **T cell receptor beta (TRB) annotation and genotype inference**. It processes **FASTA files**, infers genotypes, and OGRDB stats reports.

## Required Docker Images

Before running the pipeline, pull the necessary containers:

```bash
docker pull peresay/base:latest
docker pull williamlees/ogrdbstats:latest
```

## Running the Pipeline

### 1. Organize Input Data

Ensure your **FASTA input files** are in the correct directory. Example:

```plaintext
project_root/
├── data/
│   ├── sample1.fasta
│   ├── sample2.fasta
│   ├── sample3.fasta
```

### 2. Run the Pipeline

#### Option 1: Multiple Samples from a Directory

Process all FASTA files in a directory:

```bash
cpus=16 # Number of CPUs for nextflow to use
ncors=16 # Number of cores for IgBlast
nextflow -log logs/nextflow.log run /path/to/nextflow_demo/airrseq_annotation_pipelines/pipelines/Human_TRB_Annotation/main.nf \
  -profile docker,human_trb \
  -process.cpus ${cpus} \ 
  -w work \
  --outdir results \
  --ncors ${ncors} \ 
  --samples '/path/to/data/*.fasta'
```

#### Option 2: Multiple Samples from a CSV File

Create a CSV file with sample information:

```plaintext
sample_id,file
CI10_PRJEB28370_TRB,/path/to/data/CI10_PRJEB28370_TRB.fasta
CI11_PRJEB28370_TRB,/path/to/data/CI11_PRJEB28370_TRB.fasta
```

Then run:

```bash
cpus=16 # Number of CPUs for nextflow to use
ncors=16 # Number of cores for IgBlast
nextflow -log logs/nextflow.log run /path/to/nextflow_demo/airrseq_annotation_pipelines/pipelines/Human_TRB_Annotation/main.nf \
  -profile docker,human_trb \
  -process.cpus ${cpus} \ 
  -w work \
  --outdir results \
  --ncors ${ncors} \ 
  --input_csv samples.csv
```

#### Option 3: Single Sample

Process a single FASTA file:

```bash
cpus=16 # Total CPUs available to Nextflow
ncors=16 # Number of cores for processes
nextflow -log logs/nextflow.log run /path/to/nextflow_demo/airrseq_annotation_pipelines/pipelines/Human_TRB_Annotation/main.nf \
  -profile docker,human_trb \
  -process.cpus ${cpus} \ 
  -w work \
  --outdir results \
  --ncors ${ncors} \ 
  --fasta /path/to/data/sample.fasta
```

## Pipeline Features

### Low Depth Sample Filtering

The pipeline can automatically filter out samples with insufficient sequencing depth:

- Enable with `--check_low_depth` true (enabled by default)
- Set minimum threshold with `--check_low_depth_threshold` 100
- Samples below the threshold will be excluded from analysis with a warning message

### Genotype Inference Parameters

Customize genotype inference with:

```bash
--genotype.min_consensus_count 1  # Minimum consensus count for genotype inference
```

## Expected Output

Results are saved in the specified `--outdir` directory (default: `results/`) and include:

| File Path                                      | Description                                      |
|------------------------------------------------|--------------------------------------------------|
| `annotations/*_personal.tsv`                   | Annotated reads with personal reference set      |
| `genotype/_genotype.tsv`                       | Genotype inference results                       |
| `ogrdbstats/*_ogrdb_report.tsv`                | OGRDB statistics reports                         |
| `ogrdbstats/*_ogrdb_plots.tsv`                 | OGRDB statistics plots                           |

## Available Profiles

The pipeline includes several profiles:

- human_trb: Sets the reference directory for human TRB analysis
- docker: Runs the pipeline using Docker containers
- singularity: Runs the pipeline using Singularity containers
- singularityHPC: Runs the pipeline using Singularity on HPC with SLURM

## Debugging & Logs

To check logs:

```bash
cat logs/nextflow.log
```

For more detailed execution information:

```bash
nextflow log
```

## Resource Configuration

Adjust computational resources:

```bash
--ncors 10           # Number of cores for processes
-process.cpus 50     # Total CPUs available to Nextflow
```

## Reference Data

The pipeline requires reference data for TRB annotation:

- V_gapped.fasta
- D.fasta
- J.fasta
- human_gl.aux
- human.ndm

These are automatically set when using the human_trb profile.
