#!/bin/bash
##
## WGBS analysis (of trimmed cutadapted fastq.gz) for pe samples
##
# Izaskun Mallona
# Mon  18 2018
# GPL

TASK="cg_context_test"
DEBUG="cg_context"
WD=/home/imallona/"$TASK"
DATA="$HOME"/data
SOFT="$HOME"/soft
MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
VIRTENVS=~/virtenvs

NTHREADS=12

MAPQ_THRES=40

# taupo soft
FASTQC=/usr/local/software/FastQC/fastqc
SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
##CUTADAPT=/usr/local/bin/cutadapt
CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
BEDTOOLS="$SOFT"/bedtools/bin/bedtools

cd $WD

source $VIRTENVS/bwa-meth/bin/activate

## paired end stuff
for sample in SRR2878513_ 20151223.B-MmES_TKOD3A1c1-3_R
do
    source $VIRTENVS/bwa-meth/bin/activate

    zcat $DEBUG/"$sample"1_cutadapt_sickle.fastq.gz | head -10000 > "$sample"1_cutadapt_sickle.fastq.gz
    zcat $DEBUG/"$sample"2_cutadapt_sickle.fastq.gz | head -10000 > "$sample"2_cutadapt_sickle.fastq.gz
    
    fw="$sample"1_cutadapt_sickle.fastq.gz
    rv="$sample"2_cutadapt_sickle.fastq.gz

    ( bwameth.py --reference "$MM9" \
                 $fw $rv \
                 --threads $NTHREADS |  \
            samtools view -bS - | \
            samtools sort  - > \
                     "$sample"_bwameth_default.bam ) \
        3>&1 1>&2 2>&3 | tee "$sample"_bwameth_default.log

    deactivate

    bam="$sample"_bwameth_default.bam

    "$QUALIMAP"  bamqc \
                 -bam "$bam" \
                 -gd mm9 \
                 -outdir "$(basename $bam .bam)"_qualimap \
                 --java-mem-size=10G \
                 -nt $NTHREADS
    
    java -jar -XX:ParallelGCThreads=$NTHREADS \
         $MARKDUPLICATES INPUT=$WD/"$bam" \
         OUTPUT=$WD/"$(basename $bam .bam)""_dup_marked.bam" \
         METRICS_FILE=$WD/"$(basename $bam .bam)""_dup_marked.metrics"


    $METHYLDACKEL extract \
                  -q $MAPQ_THRES \
                  -@ $NTHREADS \
                  --cytosine_report \
                  $MM9 \
                  $bam \
                  -o $(basename $bam .bam)

    awk '
{ 
  OFS=FS="\t";
  if ($3=="+") print $1,$2-1,$2,$4,$5,$3,$7 ;
  else  print $1,$2-1,$2,$4,$5,$3,$7 ;
}
' "$(basename $bam .bam)".cytosine_report.txt  |
        "$BEDTOOLS" slop -i - \
                    -g mm9.genome \
                    -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
                    -tab \
                    -s

    # crashes due to the fact it goes outside the genome boundaries
    
    paste "$(basename $bam .bam)".cytosine_report.txt \
          "$(basename $bam .bam)"_cytosine_report_slop.fa > tmp

    ## now get odd and even lines
    awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tmp > bar
    rm -f tmp
    mv -f bar "$(basename $bam .bam)"_stranded.txt

done


