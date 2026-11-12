#!/usr/bin/env bash
# Main processing commands copied from the GM-seq pipeline used by run_spikeIn.sh.
set -euo pipefail
# Usage: bash 01_GMseq_pipeline.sh INPUT_DIR REFERENCE.fa BWA_INDEX [THREADS] [CONTROL_REFERENCE.fa CONTROL_BWA_INDEX]
work_dir=$(cd "$1" && pwd)
ref=$(realpath "$2")
bwa_index=$(realpath -m "$3")
thread=${4:-10}
control_ref=""
control_index=""
# Optional control analysis requires both CONTROL_REFERENCE.fa and CONTROL_BWA_INDEX.
if [[ $# -eq 6 ]]; then
  control_ref=$(realpath "$5")
  control_index=$(realpath -m "$6")
fi
cd "$work_dir"
mkdir -p fastqc fastqc_filt trim bam_files metrics astair logs tmp
while read -r i; do
	#QC before filtration
	fastqc -t ${thread} -o ./fastqc ./${i}.R1.fastq.gz ./${i}.R2.fastq.gz

	##filter low-quality reads and adapters
	fastp -i ./${i}.R1.fastq.gz -I ./${i}.R2.fastq.gz -o ./trim/${i}_clean.R1.fastq.gz -O ./trim/${i}_clean.R2.fastq.gz -j ./trim/${i}_fastp.json -h ./trim/${i}_fastp.html -f 10 -F 10 -q 5 -u 50 -n 14 -w ${thread} 2> ./logs/${i}_fastp.log

	##QC after filtration
	fastqc -t ${thread} -o ./fastqc_filt ./trim/${i}_clean.R1.fastq.gz ./trim/${i}_clean.R2.fastq.gz

	#BWA mapping to mm10
	bwa mem -t ${thread} -M -R "@RG\tID:${i}\tSM:${i}\tPL:ILLUMINA\tLB:${i}\tPU:${i}" ${bwa_index} ./trim/${i}_clean.R1.fastq.gz ./trim/${i}_clean.R2.fastq.gz 2> ./logs/${i}_bwa.log | samtools view -@ ${thread} -bh -o ./bam_files/${i}.bam

	#sort and index
	samtools sort -@ ${thread} -o ./bam_files/${i}.sorted.bam ./bam_files/${i}.bam
	samtools index ./bam_files/${i}.sorted.bam

	#mark/remove PCR duplicates
	picard MarkDuplicates --INPUT ./bam_files/${i}.sorted.bam --OUTPUT ./bam_files/${i}.rmdup.sorted.bam --METRICS_FILE ./metrics/${i}.markdup.txt --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT 2> ./logs/${i}_markdup.log

	samtools index -@ ${thread} ./bam_files/${i}.rmdup.sorted.bam

	##basic alignment metrics
	samtools flagstat -@ ${thread} ./bam_files/${i}.rmdup.sorted.bam > ./metrics/${i}.flagstat.txt
	samtools stats -@ ${thread} ./bam_files/${i}.rmdup.sorted.bam > ./metrics/${i}.samtools.stats.txt
	picard CollectAlignmentSummaryMetrics -R ${ref} -I ./bam_files/${i}.rmdup.sorted.bam -O ./metrics/${i}.alignment_summary.txt 2> ./logs/${i}_CollectAlignmentSummaryMetrics.log
	picard CollectInsertSizeMetrics -I ./bam_files/${i}.rmdup.sorted.bam -O ./metrics/${i}.insert_size_metrics.txt -H ./metrics/${i}.insert_size_histogram.pdf 2> ./logs/${i}_CollectInsertSizeMetrics.log

	##AsTair methylation extraction
	astair call -i ./bam_files/${i}.rmdup.sorted.bam -f ${ref} -d ./astair -t ${thread} 2> ./logs/${i}_astair.log

  # Optional spike-in mapping: original control alignment/calling commands.
  if [[ -n "$control_ref" ]]; then
    bwa mem -t ${thread} -M ${control_index} ./trim/${i}_clean.R1.fastq.gz ./trim/${i}_clean.R2.fastq.gz 2> ./logs/${i}_control_bwa.log | samtools view -@ ${thread} -bh -o ./bam_files/${i}_control.bam
    samtools sort -@ ${thread} -o ./bam_files/${i}_control.sorted.bam ./bam_files/${i}_control.bam
    samtools index ./bam_files/${i}_control.sorted.bam
    picard MarkDuplicates --INPUT ./bam_files/${i}_control.sorted.bam --OUTPUT ./bam_files/${i}_control.rmdup.sorted.bam --METRICS_FILE ./metrics/${i}_control.markdup.txt --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY LENIENT 2> ./logs/${i}_control_markdup.log
    samtools index -@ ${thread} ./bam_files/${i}_control.rmdup.sorted.bam
    astair call -i ./bam_files/${i}_control.rmdup.sorted.bam -f ${control_ref} -d ./astair -t ${thread} 2> ./logs/${i}_control_astair.log
  fi
done < filename.txt
