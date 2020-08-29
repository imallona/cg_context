#!/bin/bash
##
## new datasets single end
##
## run in taupo
##
## Filename								Description
# 20200820.A-Sample1_pool1_R1.fastq.gz 	IK_WGBS_TKO_WT
# 20200820.A-Sample2_pool1_R1.fastq.gz 	IK_WGBS_TKO_DNMT3B_9R
# 20200820.A-Sample3_pool1_R1.fastq.gz	IK_WGBS_TKO_DNMT3A2_R836K_A2
# 20200820.A-Sample4_pool1_R1.fastq.gz	IK_WGBS_TKO_DNMT3A2_R836K_D1
##
## methyldackel calls for non cpg methylation
##
## Izaskun Mallona
## 27 aug 2020
## GPL

export HOME=/home/imallona
export SRC="$HOME"/src/cg_context
export TASK="cg_context"
export WD=/home/Shared_s3it/imallona/"$TASK"/august_2020
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=30
export MAPQ_THRES=40
export MIN_DEPTH=10

export FASTQC=/usr/local/software/FastQC/fastqc
export SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
# export MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
export PICARD="$SOFT"/picard/build/libs/picard.jar
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
export FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump 

export ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
export ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"

export PROJECT=p3311
export DATA_DELIVERY=NovaSeq_20200826_NOV446_copy_o7325_DataDelivery
export CONFIG_FILE=august_2020.conf

mkdir -p $WD

cd $WD


echo "Config file"

# IL 12 and IL 18 are two biol replicates from WT DNMT3A2 in TKO
# IL 14 and 19  are two biol replicates from QC DNMT3A2 in TKO
# all samples are different clones, library prepared in paralel and sequenced in parallel
# the only difference is the DNA extraction for 12 and 14 was done in August, while 18 and 19 was done in November

cat << EOF  > "$CONFIG_FILE"
20200826.A-Sample1_pool1_R1.fastq.gz,IK_WGBS_TKO_WT
20200826.A-Sample2_pool1_R1.fastq.gz,IK_WGBS_TKO_DNMT3B_9R
20200826.A-Sample3_pool1_R1.fastq.gz,IK_WGBS_TKO_DNMT3A2_R836K_A2
20200826.A-Sample4_pool1_R1.fastq.gz,IK_WGBS_TKO_DNMT3A2_R836K_D1
EOF


echo Data retrieval

## from fabric
wget --user imallona -e robots=off --ask-password -r --no-parent -nH \
     --cut-dirs=1 \
     -A '2020*.fastq.gz' \
     --reject='*index.html*' \
     https://fgcz-gstore.uzh.ch/projects/"$PROJECT"/"$DATA_DELIVERY"


echo md5sums

cd "$PROJECT"/"$DATA_DELIVERY"

wget --user imallona -e robots=off --ask-password -r --no-parent -nH https://fgcz-sushi.uzh.ch/projects/p3311/NovaSeq_20200826_NOV446_copy_o7325_DataDelivery/md5sums.txt

md5sum -c md5sums.txt
# 20200826.A-Sample1_pool1_R1.fastq.gz: OK
# 20200826.A-Sample2_pool1_R1.fastq.gz: OK
# 20200826.A-Sample3_pool1_R1.fastq.gz: OK
# 20200826.A-Sample4_pool1_R1.fastq.gz: OK



cd $WD

echo Mapping


while IFS='' read -r line || [[ -n "$line" ]]; do
    cd "$WD"
    
    sample=$(echo $line | cut -f1 -d",")
    sample=$(basename $sample _R1.fastq.gz)
    
    r1="$(echo $sample)_R1"
    r="$r1"
    echo "$(date) Processing sample $sample"
    echo "... read 1 is $r1"

    ## process
    mkdir -p "$sample"
    cd $_
    
    ln -s "$WD"/"$PROJECT"/"$DATA_DELIVERY"/"$r1".fastq.gz
    curr="$r1"_raw
    mkdir -p $curr
    
    $FASTQC "$r".fastq.gz --outdir "$curr" \
            -t $NTHREADS &> "$curr"/"$r"_fastqc.log
    

    source $VIRTENVS/cutadapt/bin/activate
    
    cutadapt \
        -j $NTHREADS \
        -b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
        -o "$r1"_cutadapt.fastq.gz \
        "$r1".fastq.gz &> "$sample"_cutadapt.log

    deactivate

    rm -f "$r1".fastq.gz
    
    "$SICKLE" se \
              -f "$r1"_cutadapt.fastq.gz \
              -o "$r1"_cutadapt_sickle.fastq.gz \
              -t sanger \
              -g &> "$sample"_cutadapt_sickle.log


    rm -f "$r1"_cutadapt.fastq.gz 
    
    curr="$r1"_cutadapt_sickle
    mkdir -p "$curr"
    $FASTQC "$sample"_cutadapt_sickle.fastq.gz \
            --outdir "$curr" \
            -t $NTHREADS &> "$curr"/"$r"_fastqc.log
    

    source $VIRTENVS/bwa-meth/bin/activate
    
    fw="$r1"_cutadapt_sickle.fastq.gz

    bam="$sample"_bwameth_default.bam

    ( bwameth.py --reference "$MM9" \
                 "$fw" \
                 --threads $NTHREADS |  \
            samtools view --threads $NTHREADS -bS - | \
            samtools sort --threads $NTHREADS - > \
                     "$bam" ) \
        3>&1 1>&2 2>&3 | tee "$sample"_bwameth_default.log

    deactivate

    ## mapq over 40
    ## here
    samtools view -@ $NTHREADS \
             -h -b -q "$MAPQ_THRES" \
             "$bam" \
             -o "$sample"_mapq.bam
    
    mv -f "$sample"_mapq.bam "$bam"

    java -jar -XX:ParallelGCThreads=$NTHREADS \
     "$PICARD" MarkDuplicates INPUT="$bam" \
     REMOVE_DUPLICATES=TRUE \
     REMOVE_SEQUENCING_DUPLICATES=TRUE \
     OUTPUT="$(basename $bam .bam)""_dup_removed.bam" \
     METRICS_FILE="$(basename $bam .bam)""_dup_marked.metrics"

    mv -f "$bam" "$(basename $bam .bam)"with_duplicates.bam
    mv -f "$(basename $bam .bam)"_dup_removed.bam "$bam"
        
    "$QUALIMAP"  bamqc \
                 -bam "$bam" \
                 -gd mm9 \
                 -outdir "$(basename $bam .bam)"_qualimap \
                 --java-mem-size=10G \
                 -nt $NTHREADS

    cd $WD
        
done < "$CONFIG_FILE"


## tests on further stratifying the contexts
cd $WD

for bam in $(find $WD -name "*default.bam")
do
    echo $bam

    bash ~/src/cg_context/extract_motifs_frequency_from_bam.sh \
         -b $bam \
         -t $NTHREADS \
         --bedtools $BEDTOOLS \
         --methyldackel $METHYLDACKEL    
done


