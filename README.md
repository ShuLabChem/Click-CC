# Click-CC

Analysis scripts for Click-CC, which jointly profiles chromatin accessibility and drug binding in single cells.

This repository also includes scripts for bulk RNA-seq analysis and GM-seq DNA methylation processing.

## 1. Click-CC data analysis

| Script in `01_ClickCC_pipe/` | Function |
| --- | --- |
| `01_split.sh` | Separate FASTQs by modality barcode and select matching read pairs. |
| `splitbarcode0bp.pl` | Called by `01_split.sh`: exactly match the supplied barcode at a fixed R2 position and separate matching/unmatched reads. |
| `fqpick.pl` | Called by `01_split.sh`: select R1 records using read identifiers from the selected R2 FASTQ. |
| `02_pipe.sh` | Extract cell barcodes, trim and align reads, process BAMs and split by cell barcode. |
| `split.py` | Called by `02_pipe.sh`: split a BAM sorted by `CB` into one BAM per cell barcode. |
| `03_cell_annotation.R` | QC, doublet filtering, gene activity, LSI integration, clustering and marker identification from processed peak/fragment files. |

The Perl helpers separate modalities and recover read pairs; the Python helper separates aligned reads by cell. Peak calling and indexed fragment-file generation are prepared separately before the R analysis.

```bash
Rscript 01_ClickCC_pipe/03_cell_annotation.R INPUT_DIR cells.tsv ANNOTATION.rds OUTPUT_DIR
```

`cells.tsv` contains `barcode`, `sample_id` and `timepoint`. Integration uses `sample_id`. See [preprocessing commands, helper inputs/outputs and cell-analysis inputs](01_ClickCC_pipe/README.md).

## 2. RNA-seq analysis

Scripts in `02_RNA_seq/` cover read preprocessing, alignment, gene counting, differential expression analysis and comparison of drug versus labelled-drug responses.

## 3. DNA methylation processing

`03_DNA_methylation/01_GMseq_pipeline.sh` processes GM-seq reads through quality control, trimming, alignment and methylation calling. See [FASTQ format and dependencies](03_DNA_methylation/README.md).