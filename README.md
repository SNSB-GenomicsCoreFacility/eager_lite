# eager-lite

`eager-lite` is a lightweight **Nextflow pipeline** (based on nf-core/eager) designed for aligning **ancient DNA (aDNA)** sequencing data to a reference genome using **`bwa aln`** with parameters recommended for ancient samples.  
The pipeline is optimized for **large FASTQ files**, which are automatically split into smaller chunks prior to alignment to improve performance and scalability.

This workflow assumes that all sequencing libraries are **UDG-treated**.

---

## Features

- Adapter removal for single-end and paired-end reads  
- Automatic FASTQ chunking for large files  
- Alignment using `bwa aln` with ancient DNA–recommended settings  
- Library-level and sample-level BAM merging  
- Post-alignment filtering  
- PCR duplicate removal using **Picard** or **dedup**  
- Supports mixed single-end and paired-end libraries  
- Scalable and reproducible via **Nextflow**

---

## Pipeline Overview

The pipeline consists of the following major steps:

1. **Adapter Removal**  
   Sequencing adapters are removed from FASTQ files.

2. **FASTQ Splitting**  
   Large FASTQ files are split into smaller chunks to allow parallel alignment.

3. **Alignment (bwa aln)**  
   Each FASTQ chunk is aligned independently to the reference genome using `bwa aln` with settings optimized for ancient DNA.

4. **Merge Alignments by Library**  
   Chunk-level BAM files are merged per library.

5. **Merge Alignments by Sample**  
   Library-level BAM files are merged into a final sample-level BAM.

6. **Post-alignment Filtering**  
   BAM files are filtered based on mapping quality and alignment criteria.

7. **Deduplication**  
   PCR duplicates are removed using either:
   - `Picard MarkDuplicates`, or  
   - `dedup`

---

## Assumptions

- All libraries are **UDG-treated**
- Input data consists of ancient DNA FASTQ files
- Reference genome is indexed for `bwa aln`

---

## Input File Format

The pipeline requires a **comma-separated sample sheet** with the following columns:

```bash
sample,lib,udg_treated,fastq_1,fastq_2
SG24_M,lib1,true,/data/ancient/goat/SG24_M/read1_lib1.fastq.gz,/data/ancient/goat/SG24_M/read2_lib1.fastq.gz
SG24_M,lib2,true,/data/ancient/goat/SG24_M/read1_lib2.fastq.gz,none
SG2_C,lib1,true,/data/ancient/goat/SG2_C/SG2C_lib1.fastq.gz,none
```

---

## Running the Pipeline

The `eager-lite` pipeline is executed using **Nextflow** (v >25.0). At minimum, you must provide a sample sheet, a reference genome, and an output directory.

### Basic Usage

```bash
nextflow run eager_lite \
  --input samplesheet.csv \
  --fasta reference.fasta \
  --outdir results
  -profile mamba
  -qs 2
  -resume
```

---

## Output Directory Structure

The `eager-lite` pipeline generates the following directory structure inside the specified output directory:

```text
.
├── adapterremoval
├── bwa
├── dedupflagstat
├── endorspy
├── fastqc
├── filteredflagstat
├── mergelib
├── picard
├── pipeline_info
├── rawflagstat
├── samtools
├── seqkit
└── sortmergedlib

```

| Directory | Description |
|-----------|-------------|
| `adapterremoval` | Adapter-trimmed FASTQ files generated during preprocessing |
| `bwa` | BAM files produced by `bwa aln`, including chunk-level alignments |
| `mergelib` | BAM files merged at the library level from chunk-level alignments |
| `sortmergedlib` | Sorted library-level BAM files prior to sample-level merging |
| `samtools` | Sample-level merged BAM files and intermediate samtools outputs |
| `rawflagstat` | `samtools flagstat` reports for raw (unfiltered) BAM files |
| `filteredflagstat` | `samtools flagstat` reports after post-alignment filtering |
| `picard` | Deduplicated BAM files generated using Picard (when selected) |
| `dedupflagstat` | `samtools flagstat` reports for deduplicated BAM files |
| `fastqc` | FASTQC quality control reports for input and processed FASTQ files |
| `seqkit` | FASTQ chunking statistics and intermediate chunk files |
| `endorspy` | Endogenous DNA content estimation results |
| `pipeline_info` | Nextflow execution metadata, logs, reports, and software versions |

---

## References

- Apeltzer, M. **CircularMapper**
- **nf-core/eager**
- Daly, K. G. *et al.* (2021).  
  **A worldwide expansion of domesticated goats inferred from ancient DNA.**  
  *Proceedings of the National Academy of Sciences (PNAS)*, 118(20): e2100901118.  
  DOI: 10.1073/pnas.2100901118

---

## Acknowledgements

- The AdapterRemoval module used in this pipeline is copied and modified from **nf-core/eager**.
- Special thanks to the **nf-core community** for providing high-quality, open-source modules and pipelines that enable reproducible bioinformatics workflows.
