# Click-CC

<img width="1143" height="723" alt="Click-CC workflow" src="https://github.com/user-attachments/assets/eafc2364-eba0-4440-bc3a-48e14d05b68f" />

Analysis scripts for Click-CC and other related data, including bulk RNA-seq and GM-seq DNA methylation data.

## 1. Click-CC data analysis

| Script in `01_ClickCC_pipe/` | Function |
| --- | --- |
| `01_split.sh` | Separate FASTQs by modality barcode and select matching read pairs. |
| `splitbarcode0bp.pl` | Called by `01_split.sh`: exactly match the supplied barcode at a fixed R2 position and separate matching/unmatched reads. |
| `fqpick.pl` | Called by `01_split.sh`: select R1 records using read identifiers from the selected R2 FASTQ. |
| `02_pipe.sh` | Extract cell barcodes, trim and align reads, process BAMs and split by cell barcode. |
| `split.py` | Called by `02_pipe.sh`: split a BAM sorted by `CB` into one BAM per cell barcode. |
| `03_cell_annotation.R` | Generate peak-by-cell matrices with Signac, perform QC, LSI integration, clustering and marker identification. |

Run the preprocessing scripts in the corresponding paired-FASTQ directory:

```bash
bash /path/to/Click-CC/01_ClickCC_pipe/01_split.sh T7_BARCODE
bash /path/to/Click-CC/01_ClickCC_pipe/02_pipe.sh EXPECTED_CELL_NUMBER
```

After preparing peak files and indexed fragment files, run the R analysis from the repository root:

```bash
Rscript 01_ClickCC_pipe/03_cell_annotation.R INPUT_DIR cells.tsv ANNOTATION.rds OUTPUT_DIR
```

`cells.tsv` contains `barcode`, `sample_id` and `timepoint`. See [input requirements and preprocessing details](01_ClickCC_pipe/README.md).

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
