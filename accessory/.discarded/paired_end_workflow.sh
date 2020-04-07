#!/bin/bash
##
## WGBS analysis of SRA pe samples
##
# Izaskun Mallona
# Mon  18 apr 2018
# GPL

echo 'BEING REBUILT'

HOME=/home/imallona
TASK="cg_context_bulk"
CONF="$HOME"/"$TASK"/conf.sh

## this is rather stupid, just export them
cat > "$CONF" << EOL
HOME=/home/imallona
TASK="cg_context_bulk"
WD="$HOME"/"$TASK"
DATA="$HOME"/"$TASK"
SOFT="$HOME"/soft
MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
VIRTENVS=~/virtenvs

NTHREADS=18

MAPQ_THRES=40

FASTQC=/usr/local/software/FastQC/fastqc
SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
##CUTADAPT=/usr/local/bin/cutadapt
CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
BEDTOOLS="$SOFT"/bedtools/bin/bedtools
FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump 


ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"

EOL

set -a
source conf.sh
set +a

mkdir -p $WD

cd $_

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome

## fetching paired end as well
# https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP041760

## paired end stuff
for sample in SRR2878513 SRR2878520 SRR1274742 SRR1274743 SRR1274744 SRR1274745 SRR1653162
do

    echo 'uncomment to download and preprocess'
    # $FASTQDUMP -I --gzip --split-files $sample
    export sample="$sample"
    
    for r in 1 2
    do
        curr=${WD}/${sample}
        mkdir -p $curr
    
        $FASTQC ${WD}/"$sample"_"$r".fastq.gz --outdir ${curr} \
                -t $NTHREADS &> ${curr}/${sample}_fastqc.log

        source $VIRTENVS/cutadapt/bin/activate
        
        HOME=/home/imallona
        TASK="cg_context_bulk"
        CONF="$HOME"/"$TASK"/conf.sh
        source "$CONF"
        
        cutadapt \
            -j $NTHREADS \
            -b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
            -B $ILLUMINA_UNIVERSAL -B $ILLUMINA \
            -o "$WD"/"${sample}"_1_cutadapt.fastq.gz \
            -p "$WD"/"${sample}"_2_cutadapt.fastq.gz \
            "$DATA"/"$sample"_1.fastq.gz "$DATA"/"$sample"_2.fastq.gz &> "$WD"/"$sample"_cutadapt.log

        deactivate


    done

    "$SICKLE" pe \
              -f "$WD"/"$sample"_1_cutadapt.fastq.gz \
              -r "$WD"/"$sample"_2_cutadapt.fastq.gz \
              -o "$WD"/"$sample"_1_cutadapt_sickle.fastq.gz \
              -p "$WD"/"$sample"_2_cutadapt_sickle.fastq.gz \
              -t sanger \
              -s "$WD"/"$sample"_cutadapt_sickle_singles.fastq.gz \
              -g &> "$WD"/"$sample"_cutadapt_sickle.log


    rm -f "$WD"/"${sample}"_1_cutadapt.fastq.gz "$WD"/"${sample}"_2_cutadapt.fastq.gz
    
    for r in 1 2
    do
        curr="$sample"_"$r"_cutadapt_sickle
        mkdir -p "$WD"/"$curr"
        $FASTQC "$WD"/${sample}_"$r"_cutadapt_sickle.fastq.gz \
                --outdir "$WD"/"$curr" \
                -t $NTHREADS &> ${curr}/${sample}_"$r"_fastqc.log
    done    

    source $VIRTENVS/bwa-meth/bin/activate
    HOME=/home/imallona
    TASK="cg_context_bulk"
    CONF="$HOME"/"$TASK"/conf.sh
    source "$CONF"
    
    fw="$sample"_1_cutadapt_sickle.fastq.gz
    rv="$sample"_2_cutadapt_sickle.fastq.gz


    ## bwameth doesn't like reads with different read names, treats them as single end
    zcat $fw | sed 's/\.1 /\ 1 /g' | gzip -c  > "$fw"_ed.gz
    zcat $rv | sed 's/\.2 /\ 2 /g'  | gzip -c > "$rv"_ed.gz
    mv "$fw"_ed.gz  "$fw"
    mv "$rv"_ed.gz "$rv"
    
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
         REMOVE_DUPLICATES=TRUE \
         REMOVE_SEQUENCING_DUPLICATES=TRUE \
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

