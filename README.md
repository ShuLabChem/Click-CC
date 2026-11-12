# Click-CC

Analysis scripts for paired single-cell chromatin accessibility and drug binding, bulk RNA-seq, and GM-seq DNA methylation.

Update the input paths, software environments and reference files before running. Use matching genome indices and annotations; replace these files for mm10 or another genome as appropriate.

## 1. Click-CC data analysis

| Script in `01_ClickCC_pipe/` | Function |
| --- | --- |
| `01_split.sh` | Separate FASTQs by modality barcode and select matching read pairs. |
| `splitbarcode0bp.pl` | Called by `01_split.sh`: exactly match the supplied barcode at a fixed R2 position and separate matching/unmatched reads. |
| `fqpick.pl` | Called by `01_split.sh`: select R1 records using read identifiers from the selected R2 FASTQ. |
| `02_pipe.sh` | Extract cell barcodes/UMIs, trim and align reads, process BAMs and split by cell barcode. |
| `split.py` | Called by `02_pipe.sh`: split a BAM sorted by `CB` into one BAM per cell barcode. |
| `03_cell_annotation.R` | QC, doublet filtering, gene activity, LSI integration, clustering and marker identification from processed peak/fragment files. |

The Perl helpers separate modalities and recover read pairs; the Python helper separates aligned reads by cell. Peak calling and indexed fragment-file generation are prepared separately before the R analysis.

```bash
Rscript 01_ClickCC_pipe/03_cell_annotation.R INPUT_DIR cells.tsv ANNOTATION.rds OUTPUT_DIR
```

`cells.tsv` contains `barcode`, `sample_id` and `timepoint`. Integration uses `sample_id`. See [preprocessing commands, helper inputs/outputs and cell-analysis inputs](01_ClickCC_pipe/README.md).

## 2. RNA-seq analysis

Scripts in `02_RNA_seq/`:
- `01_RNA_pipe_human.sh` / `02_RNA_pipe_mouse.sh`: FASTQ QC, trimming and STAR alignment.
- `03_featureCounts.sh`: gene-level counting.
- `04_DESeq2_example.R`: differential expression analysis.
- `05_correlation_example.R`: compare log2 fold changes between drug and labelled-drug treatments.

Set the input paths and sample columns in the example scripts for your experiment.

```bash
bash 02_RNA_seq/03_featureCounts.sh ANNOTATION.gtf BAM_DIR
Rscript 02_RNA_seq/04_DESeq2_example.R DRUG_NAME
Rscript 02_RNA_seq/05_correlation_example.R DRUG_NAME
```

## 3. DNA methylation processing

`03_DNA_methylation/01_GMseq_pipeline.sh` runs the GM-seq upstream pipeline: FastQC → fastp → BWA → samtools/Picard → AsTair methylation calling.

```bash
bash 03_DNA_methylation/01_GMseq_pipeline.sh INPUT_DIR REFERENCE.fa BWA_INDEX THREADS
```

See [FASTQ format and dependencies](03_DNA_methylation/README.md).
