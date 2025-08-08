source activate tag2rose

mkdir raw 01_extract 02_align 03_bam 04_rmdup 05_merge 06_split 07_bigwig 

degree=50
cellNum=$1

cd raw
mv ../*gz ./
j=0
for i in *_R1.fastq.gz
do 
    ((j+=1))
    base=$(basename $i "_R1.fastq.gz");  
    echo $j
    umi_tools whitelist --stdin $i --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>ATAAGAGACAG{e<=2}).*'  --extract-method=regex --ed-above-threshold=correct  --error-correct-threshold=2 --set-cell-number ${cellNum} --log2stderr --knee-method=distance > ${base}_whitelist.txt --plot-prefix ${base} #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
#wait

j=0
for i in *_whitelist.txt
do
    ((j+=1))
    base=$(basename $i "_whitelist.txt")
    umi_tools extract --bc-pattern='(?P<cell_1>.{8})(?P<discard_1>ATCCACGTGCTTGAGCGCGCTGCATACTTG){e<=3}(?P<cell_2>.{6})(?P<discard_2>CCCATGATCGTCCGATCGTCGGCAGCGTCTCCACGC){e<=3}(?P<cell_3>.{6})(?P<umi_1>.{8})(?P<discard_3>ATAAGAGACAG{e<=2}).*' --stdin ${base}_R1.fastq.gz  --stdout ../01_extract/${base}_1.extract.fastq.gz  --read2-in  ${base}_R2.fastq.gz  --read2-out ../01_extract/${base}_2.extract.fastq.gz  --error-correct-cell --extract-method=regex  --whitelist ${base}_whitelist.txt #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
wait

cd ../01_extract/
#cutadapt
j=0
for i in `ls *_2.extract.fastq.gz`
do
    base=$(basename $i "_2.extract.fastq.gz")
    cutadapt -u 52 -o ${base}.trimmed.2.extract.fastq.gz ${base}_2.extract.fastq.gz
done

for i in `ls *_2.extract.fastq.gz`
do
    base=$(basename $i "_2.extract.fastq.gz")
    cutadapt -q 20 -O 10 -b CTGTCTCTTATACACATCT -b CCGAGCCCACGAGAC -B CTGTCTCTTATACACATCT -B GACGCTGCCGACGA -m 10 --max-n 0.1 --trim-n -o ${base}.cut.1.fastq.gz -p ${base}.cut.2.fastq.gz ${base}_1.extract.fastq.gz ${base}.trimmed.2.extract.fastq.gz
done

##fasta
j=0
for i in `ls *_1.extract.fastq.gz`
do
    ((j+=1))
    base=$(basename $i "_1.extract.fastq.gz")
    zcat ${base}.cut.1.fastq.gz |awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' - > ${base}_1.extract.fa #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
#wait

j=0
for i in `ls *_2.extract.fastq.gz`
do
    ((j+=1))
    base=$(basename $i "_2.extract.fastq.gz")
    zcat ${base}.cut.2.fastq.gz |awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' - > ${base}_2.extract.fa #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
#wait
cd ..

##align
j=0
for i in ./01_extract/*_2.extract.fa
do 
    ((j+=1))
    base=$(basename $i "_2.extract.fa")
    bowtie2 -p 8 -f  -N 1 --very-sensitive-local --no-unal  -x /cluster/home/jwj/reference/hg19/bowtie2/hg19 -1 ./01_extract/${base}_1.extract.fa -2 ./01_extract/${base}_2.extract.fa -S 02_align/${base}.hg19.sam 2> 02_align/${base}.hg19.align.log #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
#wait

#bam
j=0
for i in 02_align/*.sam
do
    ((j+=1))
    base=$(basename $i ".sam")
    samtools view -hbS -q 30 02_align/${base}.sam | samtools sort -T ${base} - > 03_bam/${base}.bam #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
#wait

#rmdup
j=0
for i in 03_bam/*.bam
do 
    ((j+=1))
    base=$(basename $i ".bam")
    picard MarkDuplicates -I ./03_bam/${base}.bam -O ./04_rmdup/${base}_rmdup.bam --REMOVE_DUPLICATES true -METRICS_FILE  ./04_rmdup/${base}_rmdup_picard.txt #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
#wait

#CB TB
j=0
for i in ./04_rmdup/*_rmdup.bam
do
    ((j+=1))
    base=$(basename $i "_rmdup.bam")
    samtools view ./04_rmdup/${base}_rmdup.bam | sed 's/_/\tCB:Z:/'  - | sed 's/_/\tTB:Z:/' - |awk '{print$1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$2"\t"$3}' -> ./04_rmdup/${base}_CB_TB.bam #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
wait

j=0
for i in ./04_rmdup/*_CB_TB.bam
do
    ((j+=1))
    base=$(basename $i "_CB_TB.bam")
    samtools view ./04_rmdup/${base}_rmdup.bam  -H | cat - ./04_rmdup/${base}_CB_TB.bam | samtools sort -t CB - -o ./04_rmdup/${base}_CB_TB_sorted_tags.bam
    samtools index ./04_rmdup/${base}_CB_TB_sorted_tags.bam #&
    #[ `expr $j % $degree` -eq 0 ] && wait
done
#wait

#merge
#split
for i in ./04_rmdup/*_CB_TB_sorted_tags.bam
do
    base=$(basename $i "_CB_TB_sorted_tags.bam")
    python /cluster/home/Xiaoyu/code/sc_click_epichem/sc_click_epichem_v1/split.py ./04_rmdup/${base}_CB_TB_sorted_tags.bam 06_split/ #&
done

#j=0
###bigwig
#for i in ./06_split/*.bam
#do
#    {
#    ((j+=1))
#    base=$(basename $i ".bam")
#    samtools sort -o ./06_split/${base}_sorted.bam ./06_split/${base}.bam
#    samtools index ./06_split/${base}_sorted.bam
#    num1=10000000
#    num2="$(samtools view -c ./06_split/${base}_sorted.bam  2>&1 )"
#    res=$(printf "%.5f" `echo "scale=5;$num1/$num2"|bc`)
#    bamCoverage --scaleFactor  $res -b ./06_split/${base}_sorted.bam -o ./07_bigwig/${base}.10M.bw -p 5 
#    #bamCoverage --scaleFactor  $res -b ./06_split/${base}_sorted.bam -e 300 --smoothLength 500 -o  ./07_bigwig/${base}.ext500.smo200.bw
#    } &
#    [ `expr $j % $degree` -eq 0 ] && wait
#done
#wait



