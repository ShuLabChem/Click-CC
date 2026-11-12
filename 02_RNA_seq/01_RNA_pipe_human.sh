#!/usr/bin/env bash

source activate tag2rose2

threads=$1
star_index="/cluster2/home/xiaorui/reference/hg19/STAR"

fastq_files=$(find . -name "*.fq.gz")
sample_dirs=$(dirname ${fastq_files} | sort -u)

for sample_dir in ${sample_dirs}; do
  (
    cd "${sample_dir}" 

    rename _1.fq.gz .1.fq.gz *.fq.gz
    rename _2.fq.gz .2.fq.gz *.fq.gz
    ls *.1.fq.gz | sed 's/.1.fq.gz//' > filename.txt

    mkdir -p fastqc trim star_out

    while read -r sample; do
      fastqc \
        -o fastqc \
        "./${sample}"*gz

      trim_galore \
        -q 25 \
        --phred33 \
        --length 50 \
        -e 0.1 \
        --stringency 3 \
        --paired \
        -o trim \
        "${sample}.1.fq.gz" \
        "${sample}.2.fq.gz"

      fastqc \
        -o fastqc \
        "./trim/${sample}"*gz

      STAR \
        --genomeDir "${star_index}" \
        --readFilesCommand gunzip -c \
        --readFilesIn "./trim/${sample}"*val_1* "./trim/${sample}"*val_2* \
        --runThreadN "${threads}" \
        --runMode alignReads \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic \
        --quantMode TranscriptomeSAM GeneCounts \
        --outBAMsortingThreadN 10 \
        --outFileNamePrefix "./star_out/${sample}_"
    done < filename.txt
  )
done
