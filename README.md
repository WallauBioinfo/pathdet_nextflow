# PathDet: Pathogen Detection and Classification Pipeline

**PathDet - Nextflow** is a pipeline developed in **Nextflow (DSL2)** focused on quality control, host contamination removal (*dehosting*), rapid taxonomic screening, classification via BLAST/Krona, as well as mapping and variant calling for the accurate identification of microorganisms and pathogens from NGS sequencing data (Illumina paired-end).

---

## 📋 Table of Contents

1. [Project Repository](#1-project-repository)
2. [Prerequisites and Installation](#2-prerequisites-and-installation)
3. [Data Directory Structure](#3-data-directory-structure)
4. [Input File Format (CSV)](#4-input-file-format-csv)
5. [Configuration Parameters](#5-configuration-parameters)
6. [Resource Management](#6-resource-management)
7. [How to Run](#7-how-to-run)
8. [Workflows](#8-workflows)
9. [Output Structure](#9-output-structure)

---

## 1. Project Repository

To clone and get started with the project in your local machine or HPC environment:

```bash
git clone https://github.com/WallauBioinfo/pathdet_nextflow.git
cd pathdet_nextflow
```

---

## 2. Prerequisites and Installation

Create the Conda environment from `env/env.yml`:

```bash
micromamba env create -f env/env.yml
micromamba activate pathdet_env
```

Two tools are not installed by the package manager and must be set up once, with the environment active.

### SparK

The coverage plots generated per contig (`ALIGN_AND_VARIANT_CALLING`) use [SparK](https://github.com/harbourlab/SparK) (not Apache Spark). The pipeline expects it at `$CONDA_PREFIX/opt/git/SparK/SparK.py`, configurable with `--spark_script`.

```bash
git clone https://github.com/harbourlab/SparK.git $CONDA_PREFIX/opt/git/SparK
```

### KronaTools taxonomy

`MERGE_AND_KRONA` needs the KronaTools taxonomy and the accession-to-taxid table. Without `all.accession2taxid.sorted`, the step fails with `Sorted accession to taxID list not found`.

```bash
bash $CONDA_PREFIX/opt/krona/updateTaxonomy.sh
bash $CONDA_PREFIX/opt/krona/updateAccessions.sh
```

Both scripts download data from NCBI (several GB for the accessions), so run them on a node with internet access. The files are created in `$CONDA_PREFIX/opt/krona/taxonomy/`.

### Hardware

The size of the databases defines the requirements:

| Resource | Reference |
| :--- | :--- |
| `blast/nt` | about 719 GB on disk; memory-mapped by BLASTN, so a large page cache speeds it up |
| `kraken2/k2_pluspf_20241228` | about 92 GB on disk; `RAPID_ANALYSIS` loads it entirely into RAM (about 90 GB) |

---

## 3. Data Directory Structure

> **Important:** The data directory (`params.datadir`) containing the taxonomic databases and reference indices **must be provided by the user**. The Nextflow pipeline assumes that this directory is already structured and ready before execution.

The pipeline looks for the data directory provided via the `--datadir` parameter. The expected structure inside this directory is:

```text
db/
├── blast/
│   ├── nt/                            # BLASTN (nt) database for identification
│   └── t2t/                           # BLASTN database for additional dehosting
├── bowtie2/
│   └── GRCh38_noalt_as/               # Host Bowtie2 index
├── hisat2/
│   └── grch38_tran/                   # Host HISAT2 index
├── kraken2/
│   ├── human/                         # Kraken2 index for human dehosting
│   └── k2_pluspf_20241228/            # Kraken2 database for rapid screening (PathDB)
├── ncbi_genome/
│   ├── all/                           # Full genome database files
│   ├── reference/                     # Reference genome database files
│   └── ref_list.txt                   # Reference genome list file
└── taxonkit/
    └── taxdump/                       # NCBI taxonomy files (nodes.dmp, names.dmp)
```

To prepare the database, follow section 2 of the following document: https://github.com/khoriba/PATHDET_long-read/blob/main/README.md

---

## 4. Input File Format (CSV)

The pipeline accepts a CSV file specifying the sample IDs and paths to the FASTQ files (`paired-end`).

Example (`assets/reads/reads_test.csv`):

```csv
sample_id,reads_r1,reads_r2
sample_01,path/to/sample_01_R1.fastq.gz,path/to/sample_01_R2.fastq.gz
sample_02,path/to/sample_02_R1.fastq.gz,path/to/sample_02_R2.fastq.gz
```

---

## 5. Configuration Parameters

Parameters are set in `nextflow.config` or on the command line (`--parameter value`).

### Files and Directories
| Parameter | Description | Default Value |
| :--- | :--- | :--- |
| `--input_csv` | Path to the input CSV file | `${launchDir}/assets/reads/reads_test.csv` |
| `--out` | Output directory for results | `${launchDir}/results` |
| `--datadir` | Path to the directory containing databases | `${launchDir}/assets/db` |
| `--cleanup_work_on_success` | Deletes the `work/` directory when the run succeeds. Set to `false` to use `-resume` | `true` |

### Quality Control (Fastp and BBDuk)
| Parameter | Description | Default Value |
| :--- | :--- | :--- |
| `--fastp_qualified_phred` | Minimum Phred quality score accepted by Fastp | `30` |
| `--fastp_n_base_limit` | Maximum number of uncalled bases (N) | `20` |
| `--fastp_length_required` | Minimum read length required | `80` |
| `--fastp_trim_front1` | Bases trimmed from the 5' end of read 1 | `15` |
| `--fastp_trim_tail1` | Bases trimmed from the 3' end of read 1 | `1` |
| `--fastp_trim_tail2` | Bases trimmed from the 3' end of read 2 | `1` |
| `--bbduk_entropy` | Entropy threshold for low-complexity sequence filtering | `0.5` |
| `--bbduk_entropywindow` | BBDuk entropy window size | `10` |
| `--bbduk_entropyk` | BBDuk entropy k-mer size | `5` |

### BLAST
| Parameter | Description | Default Value |
| :--- | :--- | :--- |
| `--split_size` | Number of FASTA **lines** per chunk (2 lines per read, so `50000` is 25,000 reads) | `50000` |
| `--blast_batch_chunks` | Maximum number of chunks, from all samples, searched together in one BLASTN run | `50` |

### Databases and Tools
Derived from `--datadir`; each one can be overridden individually.

| Parameter | Default Value |
| :--- | :--- |
| `--kraken2_grch38` | `<datadir>/kraken2/human` |
| `--bowtie2_host` | `<datadir>/bowtie2/GRCh38_noalt_as/GRCh38_noalt_as` |
| `--hisat2_grch38r` | `<datadir>/hisat2/grch38_tran/genome_tran` |
| `--blast_host2` | `<datadir>/blast/t2t/GCF_009914755.1_T2T-CHM13v2.0_genomic` |
| `--pathdb1_db` | `<datadir>/kraken2/k2_pluspf_20241228` |
| `--blast_nt` | `<datadir>/blast/nt/nt` |
| `--ncbi_genome_dir` | `<datadir>/ncbi_genome` |
| `--ncbi_lineage_tsv` | `<datadir>/ncbi_genome/reference/refseq_lineage_taxid_fasta.tsv` |
| `--spark_script` | `$CONDA_PREFIX/opt/git/SparK/SparK.py` |

### Resources
| Parameter | Description | Default Value |
| :--- | :--- | :--- |
| `--max_cpus` | Maximum CPUs the pipeline may use | Detected (see [Resource Management](#6-resource-management)) |
| `--max_memory` | Maximum memory in GB | 90% of the machine memory |
| `--threads_<PROCESS>` | CPUs/threads of one process | See the table in section 6 |
| `--mem_<PROCESS>` | Memory in GB of one process | See the table in section 6 |

---

## 6. Resource Management

No thread count is fixed inside the modules: every tool receives the CPUs of its process, so changing a process changes all the tools it runs.

### Limits

`--max_cpus` and `--max_memory` cap everything the pipeline requests, and the local executor never schedules more than these values, so heavy tasks wait instead of oversubscribing the machine. When not provided:

- CPUs are read from `PBS_NP`, `PBS_NUM_PPN`, `NCPUS` or `SLURM_CPUS_ON_NODE`, falling back to the CPUs of the machine.
- Memory is 90% of the machine memory.

Set them explicitly when the scheduler allocation is smaller than the node:

```bash
nextflow run main.nf --max_cpus 20 --max_memory 170
```

### Defaults per process

Values are capped by `--max_cpus` and `--max_memory`.

| Process | Threads | Memory (GB) |
| :--- | :--- | :--- |
| `QUALITY_CONTROL` | 8 | 8 |
| `DEHOST_KRAKEN2` | 8 | database size + 4 (minimum 8) |
| `DEHOST_BOWTIE2` | 8 | 8 |
| `DEHOST_HISAT2` | 8 | 12 |
| `DEHOST_BLASTN` | 4 | 4 |
| `RAPID_ANALYSIS` | 4 | database size + 8 (minimum 16) |
| `SPLIT_FASTA` | 1 | 4 |
| `RUN_BLASTN` | half of `--max_cpus` | 60 |
| `MERGE_AND_KRONA` | 2 | 8 |
| `FETCH_REFERENCE_GENOMES` | 1 | 4 |
| `ALIGN_AND_VARIANT_CALLING` | 4 | 8 |

A task that fails with exit status 137, 140, 143 or 247 is retried up to 2 times with the memory multiplied by the attempt number.

### Changing one process

Use the process name on the command line:

```bash
nextflow run main.nf --threads_RUN_BLASTN 8 --mem_RAPID_ANALYSIS 120
```

---

## 7. How to Run

Run the complete pipeline:

```bash
nextflow run main.nf \
  --input_csv "/path/to/your_samples.csv" \
  --datadir "/path/to/your/databases/db"
```

To resume an interrupted execution while reusing completed steps, keep the `work/` directory and add `-resume`:

```bash
nextflow run main.nf \
  --input_csv "/path/to/your_samples.csv" \
  --datadir "/path/to/your/databases/db" \
  --cleanup_work_on_success false \
  -resume
```

To save execution reports:

```bash
nextflow run main.nf \
  --input_csv "/path/to/your_samples.csv" \
  --datadir "/path/to/your/databases/db" \
  -with-report report.html \
  -with-timeline timeline.html \
  -with-trace trace.txt
```

---

## 8. Workflows

1. **Quality Control (QC)**: Quality filtering, adapter trimming via `fastp`, and low-entropy sequence removal via `bbduk`.
2. **Dehosting**: Sequential filtering of human/host reads using `Kraken2`, `Bowtie2`, `HISAT2`, and `BLASTN`.
3. **Rapid Analysis**: Fast k-mer based screening (`Kraken2`), taxonomic diversity reports, and an interactive Krona HTML plot.
4. **BLAST & Krona**: Splitting filtered FASTA files into chunks (`split_fasta.py`), searching the chunks of all samples together in a single `BLASTN` run against the nucleotide database, splitting the result back per sample (`blast_batch.py`), aggregating hits, and building detailed reports (`process_blast_report.py`).
5. **Mapping and Variants**: Automatic search for reference genomes based on top taxonomic hits (`fetch_genomes.py`), alignment via `minimap2`, variant calling via `bcftools`, and coverage plots using `bamCoverage` and `SparK.py`.

A `BLASTN` search against `nt` costs about one pass over the database regardless of the number of reads, so searching all samples together takes roughly the time of a single sample. Each sample still receives exactly the result it would get when searched alone.

---

## 9. Output Structure

Results are written to `--out`, one directory per sample:

```text
results/
└── <sample_id>/
    ├── <sample_id>_qc/       # Quality control
    ├── <sample_id>_hgs/      # Dehosting
    ├── <sample_id>_rapid/    # Rapid Kraken2 screening reports
    ├── <sample_id>_blast/    # BLAST results and Krona plots
    ├── <sample_id>_tbl/      # BLAST summary tables and reports
    └── <sample_id>_map/      # Reference genomes, alignments, variants and coverage
```
