#!/bin/bash
T7=$1

mkdir raw
mkdir 02_${T7}

pair_reads(){
base=$(basename $1 "_R2.fastq.gz")

gunzip ${base}_R1.fastq.gz &
gunzip ${base}_R2.fastq.gz &
wait

#sed '/^@/s/\/1$//g' ${base}_R1.fastq > raw/${base}_R1.fastq &
#sed '/^@/s/\/2$//g' ${base}_R2.fastq > raw/${base}_R2.fastq &
awk 'NR%4==1 {gsub("/1$", ""); print} NR%4!=1 {print}' ${base}_R1.fastq > raw/${base}_R1.fastq &
awk 'NR%4==1 {gsub("/2$", ""); print} NR%4!=1 {print}' ${base}_R2.fastq > raw/${base}_R2.fastq &
wait

gzip raw/${base}_R1.fastq &
gzip raw/${base}_R2.fastq &
wait
}
export -f pair_reads
ls *_R2.fastq.gz | parallel -j20 pair_reads

#rm -rf *fastq

cd ./raw

split_modality() {
    T7=$1
    base=$(basename $2 "_R2.fastq.gz")
    # short T7-21 TCCAAG and prepare splitbarcode0bp.pl file 
    perl /cluster/home/Xiaoyu/code/sc_click_epichem/sc_click_epichem_v1/splitbarcode0bp.pl ../raw/${base}_R2.fastq.gz 27 ${T7} ../02_${T7}/${base}_R2.fastq.gz ${base}_${T7}_noid.2.fastq.gz

}
export -f split_modality
ls *_R2.fastq.gz | parallel -j20 split_modality ${T7}


cd ../
pick_read1() {
    T7=$1
    i=$2
    base=$(basename $2 "_R2.fastq.gz")
    perl /cluster/home/Xiaoyu/code/sc_click_epichem/sc_click_epichem_v1/fqpick.pl $i ./raw/${base}_R1.fastq.gz ./02_${T7}/${base}_R1.fastq.gz
}
export -f pick_read1
ls ./02_${T7}/*_R2.fastq.gz | parallel -j30 pick_read1 ${T7}

