# PathDet: Pathogen Detection and Classification Pipeline

**PathDet - Nextflow** is a pipeline developed in **Nextflow (DSL2)** focused on quality control, host contamination removal (*dehosting*), rapid taxonomic screening, classification via BLAST/Krona, as well as mapping and variant calling for the accurate identification of microorganisms and pathogens from NGS sequencing data (Illumina paired-end).

---

## 📋 Table of Contents

1. [Project Repository](#1-project-repository)
2. [Prerequisites and Installation](#2-prerequisites-and-installation)
3. [Data Directory Structure](#3-data-directory-structure)
4. [Input File Format (CSV)](#4-input-file-format-csv)
5. [Configuration Parameters](#5-configuration-parameters)
6. [How to Run](#6-how-to-run)
7. [Workflows](#7-workflows)

---

## 1. Project Repository

To clone and get started with the project in your local machine or HPC environment:

```bash
git clone https://github.com/WallauBioinfo/pathdet_nextflow.git
cd pathdet_nextflow
```

---

## 2. Prerequisites and Installation

You can create the Conda environment using the `environment.yml` file:

```bash
micromamba env create -f env/env.yml
micromamba activate pathdet_envg
```

### Apache Spark & SparK.py
The coverage plot generation script per contig (`ALIGN_AND_VARIANT_CALLING`) requires **Spark** to be installed and accessible within the environment (typically mapped at `$CONDA_PREFIX/opt/git/SparK/SparK.py`).

```bash
git clone https://github.com/harbourlab/SparK.git $CONDA_PREFIX/opt/git/SparK
bash $CONDA_PREFIX/opt/krona/updateTaxonomy.sh
```

---

## 3. Data Directory Structure

> **Important:** The data directory (`params.datadir`) containing the taxonomic databases and reference indices **must be provided by the user**. The Nextflow pipeline assumes that this directory is already structured and ready before execution.

the pipeline looks for the data directory provided via the `--datadir` parameter. The expected structure inside this directory is:

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

Below are the main configurable parameters in `nextflow.config` or via command line (`--parameter value`):

### Files and Directories
| Parameter | Description | Default Value |
| :--- | :--- | :--- |
| `--input_csv` | Path to the input CSV file | `${launchDir}/assets/reads/reads_test.csv` |
| `--out` | Output directory for results | `${launchDir}/results` |
| `--datadir` | Path to the directory containing databases | None |
| `--cleanup_work_on_success` | Cleans up the `work/` directory upon success | `true` |

### Quality Control (QC: Fastp and BBDuk)
| Parameter | Description | Default Value |
| :--- | :--- | :--- |
| `--fastp_qualified_phred` | Minimum Phred quality score accepted by Fastp | `30` |
| `--fastp_n_base_limit` | Maximum limit of uncalled bases (N) | `20` |
| `--fastp_length_required` | Minimum read length required | `80` |
| `--bbduk_entropy` | Entropy threshold for low-complexity sequence filtering | `0.5` |

### BLAST Processing and Mapping
| Parameter | Description | Default Value |
| :--- | :--- | :--- |
| `--split_size` | Number of reads per FASTA batch sent to BLAST | `50000` |
| `--blast_threads` | Threads allocated for BLASTN executions | `20` |
| `--map_threads` | Threads allocated for Minimap2 alignment | `16` |

---

## 6. How to Run

To run the complete pipeline enabling the Conda profile:

```bash
nextflow run main.nf \
  --input_csv "/path/to/your_samples.csv" \
  --datadir "/path/to/your/databases/db" \
```

To resume an interrupted execution while reusing completed steps, use the `-resume` flag:

```bash
nextflow run main.nf \
  --input_csv "/path/to/your_samples.csv" \
  --datadir "/path/to/your/databases/db" \
  -resume
```

---

## 7. Workflows

1. **Quality Control (QC)**: Quality filtering, adapter trimming via `fastp`, and low-entropy sequence removal via `bbduk`.
2. **Dehosting**: Sequential filtering of human/host reads using `Kraken2`, `Bowtie2`, `HISAT2`, and `BLASTN`.
3. **Rapid Analysis**: Fast k-mer based screening (`Kraken2`), taxonomic diversity reports, and an interactive Krona HTML plot.
4. **BLAST & Krona**: Splitting filtered FASTA files into batches (`split_fasta.py`), aligning against the nucleotide database (`BLASTN`), aggregating hits, and building detailed reports (`process_blast_report.py`).
5. **Mapping and Variants**: Automatic search for reference genomes based on top taxonomic hits (`fetch_genomes.py`), alignment via `minimap2`, variant calling via `bcftools`, and coverage plots using `bamCoverage` and `SparK.py`.