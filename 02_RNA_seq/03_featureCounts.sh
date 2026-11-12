#!/usr/bin/env bash

source activate tag2rose2

# Usage: bash 03_featureCounts.sh <annotation.gtf> <bam_dir> [bam_dir ...]

gtf=$1
shift
threads=${THREADS:-16}

mapfile -t bam_files < <(
  find "$@" -type f -name "*Aligned.sortedByCoord.out.bam" | sort
)

featureCounts -T "${threads}" -a "${gtf}" -o RNA_raw_count.txt -p -B -C -t exon -g gene_id "${bam_files[@]}"

sample_count=${#bam_files[@]}
last_count_column=$((sample_count + 6))

awk 'NR > 1' RNA_raw_count.txt |
  cut -f "1,6-${last_count_column}" \
  > RNA_count_and_length.tmp

sed -n '2p' RNA_raw_count.txt |
  cut -f 1,6- |
  awk -F '\t' '
    BEGIN {
      OFS = "\t"
    }
    {
      printf "%s\t%s", $1, $2
      for (i = 3; i <= NF; i++) {
        n = split($i, path, "/")
        sample = path[n]
        sub("_Aligned.sortedByCoord.out.bam", "", sample)
        printf "\t%s", sample
      }
      printf "\n"
    }
  ' > RNA_header.tmp

sed -n '2,$p' RNA_count_and_length.tmp > RNA_counts.tmp
cat RNA_header.tmp RNA_counts.tmp > RNA_count_and_length.txt

rm RNA_header.tmp RNA_counts.tmp RNA_count_and_length.tmp
