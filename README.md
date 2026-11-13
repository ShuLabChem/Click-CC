# Click-CC

<img width="1143" height="723" alt="Click-CC workflow" src="https://github.com/user-attachments/assets/eafc2364-eba0-4440-bc3a-48e14d05b68f" />

Analysis scripts for Click-CC and other related data, including bulk RNA-seq and GM-seq DNA methylation data.

## 1. Click-CC data analysis

| Script in `01_ClickCC_pipe/` | Function |
| --- | --- |
| `01_split.sh` | Standardize paired-read identifiers, separate FASTQs by modality barcode and select matching read pairs. |
| `splitbarcode0bp.pl` | Called by `01_split.sh`: exactly match the supplied barcode at a fixed R2 position and write matching/unmatched reads to separate gzipped FASTQs. |
| `fqpick.pl` | Called by `01_split.sh`: select R1 records using read identifiers from the selected R2 FASTQ. |
| `02_pipe.sh` | Extract cell barcodes/UMIs, trim reads, align with Bowtie2, process BAMs and split by cell barcode. |
| `split.py` | Called by `02_pipe.sh`: split a BAM sorted by `CB` into one BAM per cell barcode. |
| `03_cell_annotation.R` | Build ATAC/drug matrices, perform QC, doublet filtering, gene activity, LSI integration, clustering and marker ranking. |

### 1.1. Separate modalities and select read pairs

Run `01_split.sh` in the directory containing `PREFIX_R1.fastq.gz` and `PREFIX_R2.fastq.gz`. Update the helper-script paths in the shell script to their installed locations.

```bash
bash /path/to/Click-CC/01_ClickCC_pipe/01_split.sh T7_BARCODE
```

`T7_BARCODE` is the modality-identifying sequence. The script removes terminal `/1` and `/2` mate suffixes from FASTQ headers, writes the normalized files to `raw/`, and calls the two Perl helpers:

- **`splitbarcode0bp.pl`** compares R2 with the supplied barcode starting at zero-based position 27 (the 28th base), allowing zero mismatches. It writes complete reads without trimming the barcode. Selected R2 reads go to `02_<T7_BARCODE>/`; unmatched R2 reads go to the `*_noid.2.fastq.gz` file under `raw/`.
- **`fqpick.pl`** collects the selected R2 read identifiers using the first whitespace-delimited field of each FASTQ header, then retains matching R1 records in their original R1 order. Selected R1 reads go to `02_<T7_BARCODE>/`.

Individual helper usage, with script paths adjusted to their installed locations:

```bash
perl splitbarcode0bp.pl INPUT_R2.fastq.gz 27 T7_BARCODE SELECTED_R2.fastq.gz UNMATCHED_R2.fastq.gz
perl fqpick.pl SELECTED_R2.fastq.gz INPUT_R1.fastq.gz SELECTED_R1.fastq.gz
```

### 1.2. Extract cell barcodes and process alignments

Run `02_pipe.sh` in the directory containing the paired FASTQs selected for one modality:

```bash
bash /path/to/Click-CC/01_ClickCC_pipe/02_pipe.sh EXPECTED_CELL_NUMBER
```

The script uses UMI-tools to build a barcode whitelist and extract three cell-barcode segments (8 + 6 + 6 bases) and an 8-base UMI from R1. Cutadapt trims the reads, Bowtie2 aligns them, and samtools/Picard process the alignments. The script adds `CB` and `TB` tags from the read names, groups alignments by `CB`, and calls `split.py` to write per-cell BAMs to `06_split/`.

`split.py` requires reads with the same cell barcode to be contiguous. Its input is a BAM sorted by `CB`; the output directory must exist and end with `/`. To run it separately, use:

```bash
samtools sort -t CB INPUT.bam -o CB_sorted.bam
mkdir -p per_cell
python /path/to/Click-CC/01_ClickCC_pipe/split.py CB_sorted.bam per_cell/
```

Preprocessing requires GNU parallel, Perl, Python/pysam, UMI-tools, cutadapt, Bowtie2, samtools and Picard. Update the software environment, Bowtie2 index and helper-script paths in the shell scripts. The supplied alignment example uses hg19; use the matching mm10 index for mouse data.

Peak calling and indexed fragment-file generation are separate preparation steps before the following R analysis.

### 1.3. Cluster cells and identify markers

The input directory contains:

```text
bulk/04_macs2/*.ATAC_peaks.narrowPeak
bulk/04_macs2/*.drug_peaks.narrowPeak
sinto/06_fragTab/ATAC_merge.sort.tsv.gz
sinto/06_fragTab/drug_merge.sort.tsv.gz
```

Provide one peak file per matching pattern and the tabix index for each fragment file. `drug` denotes the drug-binding channel. `ANNOTATION.rds` is a Signac-compatible gene annotation for the same genome. Run from the repository root:

```bash
Rscript 01_ClickCC_pipe/03_cell_annotation.R INPUT_DIR cells.tsv ANNOTATION.rds OUTPUT_DIR [DATASET_ID]
```

`cells.tsv` contains `barcode`, `sample_id` and `timepoint`; see [cell_metadata.example.tsv](01_ClickCC_pipe/cell_metadata.example.tsv). Barcodes must match the fragment files. Integration uses `sample_id`, and QC/display use `timepoint`. The retained QC settings cover Day1 and Day12. The optional `DATASET_ID` names the input library and defaults to `dataset`.

R dependencies include Seurat, Signac, rtracklayer, tidyverse, scDblFinder, patchwork, qs, data.table and COSG. Main outputs are `clusters.pdf`, `cluster_marker_genes.csv` and `Rds/obj_clustered.qs`, for manual cell-type annotation.

## 2. RNA-seq analysis

Scripts in `02_RNA_seq/` cover read preprocessing, alignment, gene counting, differential expression and drug versus labelled-drug comparisons. After preprocessing and alignment, run from the repository root, with paths and sample groups configured in the R examples:

```bash
bash 02_RNA_seq/03_featureCounts.sh ANNOTATION.gtf BAM_DIR
Rscript 02_RNA_seq/04_DESeq2_example.R DRUG_NAME
Rscript 02_RNA_seq/05_correlation_example.R DRUG_NAME
```

## 3. DNA methylation processing

The GM-seq pipeline performs read preprocessing, alignment and methylation calling. Prepare paired FASTQs and a sample-prefix list (`filename.txt`) in `INPUT_DIR`, then run from the repository root:

```bash
bash 03_DNA_methylation/01_GMseq_pipeline.sh INPUT_DIR REFERENCE.fa BWA_INDEX THREADS
```

See [FASTQ format and dependencies](03_DNA_methylation/README.md).
