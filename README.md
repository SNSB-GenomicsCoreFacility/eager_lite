# Ancient mtDNA Assembly Pipeline (Circularized Reference Approach)

## Overview

This **Nextflow** pipeline performs **mitochondrial DNA (mtDNA) assembly** from sequencing data of ancient samples using a **circularized reference approach**.  
It is designed to generate consensus sequences of mtDNA by circularizing the reference genome, aligning reads, and assembling using ANGSD-based consensus methods optimized for ancient DNA.

## Workflow Summary

1. **Adapter Removal**
   - Raw sequencing reads are cleaned using the **AdapterRemoval** module.
   - This module is **directly copied and modified** from the [nf-core/eager](https://github.com/nf-core/eager) pipeline to ensure compatibility with ancient DNA data and consistency in preprocessing.

2. **Reference Circularization**
   - The pipeline creates a **circularized reference** genome based on a user-provided number of base pairs (`--circle_nbp`).
   - This is achieved by appending and prepending the specified number of bases to simulate the circular nature of mtDNA.

3. **Read Alignment**
   - The circularized reference is indexed and used to align the sequencing reads with **BWA aln**, optimized for ancient DNA reads.

4. **Realignment of Decircularized Reference**
   - The aligned reads are processed using the **`realignsamfile`** tool from [CircularMapper](https://github.com/apeltzer/CircularMapper).
   - This step generates a BAM file corresponding to the **de-circularized** reference.

5. **Consensus Sequence Generation**
   - **ANGSD** is used to generate consensus sequences from the aligned BAM files.
   - The ANGSD settings follow the recommendations from the study:  
     [*A worldwide expansion of domesticated goats inferred from ancient DNA*](https://www.pnas.org/doi/10.1073/pnas.2100901118).

## Dependencies

Ensure the following tools are available in your environment or installed via your Nextflow profile:

- [Nextflow](https://www.nextflow.io/) v 25.04
- conda/mamba

## Input Requirements

- **FASTQ file list:**  
  A CSV file containing paths to the paired-end FASTQ files. Example:

```bash
sample,lib,udg_treated,fastq1,fastq2
sample1,lib1,true,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample1,lib2,false,/path/to/sample1_R1.fastq.gz,none
```
## Example Command

Run the pipeline as follows (when running on interactive node such as cm4_inter in LRZ linux cluster):

```bash
nextflow run assemble_ancient_mtDNA/ \
--input sample_fastq.2.csv \
--fasta ../references/capra_hircus/genome.fa \
--bwa_idx ../references/capra_hircus/bwa_idx/
--outdir $SCRATCH/goat/ancient_samples/screening_run \
-profile lrz_serial_std,mamba \
-resume \
-w $SCRATCH
```

- **important notes** when running on HPC with multiple clusters such as LRZ linux cluster the default partition is "serial_std". Therefore, first, type the following in your environment:

```bash
export SLURM_CLUSTERS=serial
```
then, run the above-mentioned command.


## References
- Apeltzer, M. [CircularMapper](https://github.com/apeltzer/CircularMapper)
- [nf-core/eager](https://nf-co.re/eager/2.5.2/)
- Daly, K. G. et al. (2021). A worldwide expansion of domesticated goats inferred from ancient DNA.
Proceedings of the National Academy of Sciences (PNAS), 118(20): e2100901118.
DOI: 10.1073/pnas.2100901118

## Acknowledgements
- The AdapterRemoval module used in this pipeline is copied and modified from the nf-core/eager
- Special thanks to the nf-core community for providing high-quality, open-source modules and pipelines that enable reproducible bioinformatics workflows.
