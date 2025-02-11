# Nextflow AIRR-seq Annotation Pipelines Repository

## Overview

This repository contains multiple Nextflow-based pipelines for adaptive immune receptor sequences (AIRR-seq) annotation and analysis. Each pipeline is structured to process input sequencing data, perform annotation, and generate relevant reports.

---

## Available Pipelines

### **1. TRB Annotation Pipeline**

- **Description**: Annotates T cell receptor beta (TRB) sequences, infers genotypes, and generates OGRDB stats reports.
- **Path**: `airrseq_annotation_pipelines/pipelines/Human_TRB_Annotation/`
- **Usage**: See the [TRB Annotation Pipeline README](https://github.com/ayeletperes/airrseq_annotation_pipelines/blob/main/pipelines/Human_TRB_Annotation/README.md) for details.

---

Part of the module code was adapted from [the flairr dsl2 script](https://github.com/williamdlees/flairr_dsl2) from written by William Lees

## Installation

The pipeline will run on Linux, or Windows Subsystem for Linux (WSL).

Required installations:

### **1. Install Nextflow**

[Nextflow](https://www.nextflow.io/) requires Java 8 or later. 

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow /usr/local/bin/
```

or

```bash
conda install -c bioconda nextflow
```

### **2. Install Docker**

This pipeline runs in a [**Docker container**](https://www.docker.com/).

```bash
sudo apt update && sudo apt install docker.io
sudo systemctl start docker
sudo systemctl enable docker
```

---

## Required Docker Images

Before running the pipelines, pull the necessary containers:

```bash
docker pull peresay/base:latest
docker pull williamlees/ogrdbstats:latest
```

---

## Running a Pipeline

Each pipeline has its own **README.md** inside its respective directory. Navigate to the specific pipeline and follow the instructions.

Example for the **Human TRB Annotation Pipeline**:

```bash
cd airrseq_annotation_pipelines/pipelines/Human_TRB_Annotation/
cat README.md
```

---
