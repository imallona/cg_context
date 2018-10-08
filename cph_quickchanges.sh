#!/bin/bash
##
## methyldackel calls for non cpg methylation
## outputs are processed by cph_calls.Rmd
##
# Izaskun Mallona
# aug the 22nd
# GPL

export HOME=/home/imallona
export TASK="cg_context"
export WD="$HOME"/mnt/nfs/"$TASK"/quickchanges
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=24
export MIN_DEPTH=10

export FASTQC=/usr/local/software/FastQC/fastqc
export SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
export FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump 

export ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
export ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"

mkdir -p $WD

cd $WD


# cd $WD

## no idea which one is each! @ todo ask 
cat << EOF  > quickchanges.conf
20180913.A-WGBS_3,20180913.A-WGBS_3_R1,20180913.A-WGBS_3_R2,bwa_hiseq4k_pe
20180913.A-WGBS_4,20180913.A-WGBS_4_R1,20180913.A-WGBS_4_R2,bwa_hiseq4k_pe
EOF

echo data retrieval

## from fabric
wget --user imallona -e robots=off --ask-password -r --no-parent -nH \
     --cut-dirs=2 \
     -A '*WGBS*gz' \
     --reject='index.html*' \
     https://fgcz-gstore.uzh.ch/projects/p2046/HiSeq4000_20180913_RUN481_o4758_DataDelivery/


echo next

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome



while IFS='' read -r line || [[ -n "$line" ]]; do
    cd "$WD"
    
    sample=$(echo $line | cut -f1 -d",")
    r1=$(echo $line | cut -f2 -d",")
    r2=$(echo $line | cut -f3 -d",")

    echo "$(date) Processing sample $sample"
    echo "... read 1 is $r1"
    echo "... read 2 is $r2"

    ## process

    for r in "$r1" "$r2"
    do
        ln -s "HiSeq4000_20180913_RUN481_o4758_DataDelivery"/"$r"
        curr="$r"_raw
        mkdir -p $curr
        
        $FASTQC "$r".fastq.gz --outdir "$curr" \
                -t $NTHREADS &> "$curr"/"$r"_fastqc.log
    done

    source $VIRTENVS/cutadapt/bin/activate
    
    cutadapt \
        -j $NTHREADS \
        -b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
        -B $ILLUMINA_UNIVERSAL -B $ILLUMINA \
        -o "$r1"_cutadapt.fastq.gz \
        -p "$r2"_cutadapt.fastq.gz \
        "$r1".fastq.gz "$r2".fastq.gz &> "$sample"_cutadapt.log

    deactivate

    rm -f "$r1".fastq.gz rm -f "$r2".fastq.gz
    
    "$SICKLE" pe \
              -f "$r1"_cutadapt.fastq.gz \
              -r "$r2"_cutadapt.fastq.gz \
              -o "$r1"_cutadapt_sickle.fastq.gz \
              -p "$r2"_cutadapt_sickle.fastq.gz \
              -t sanger \
              -s "$sample"_cutadapt_sickle_singles.fastq.gz \
              -g &> "$sample"_cutadapt_sickle.log


    rm -f "$r1"_cutadapt.fastq.gz "$r2"_cutadapt.fastq.gz
    
    for r in "$r1" "$r2"
    do
        curr="$r"_cutadapt_sickle
        mkdir -p "$curr"
        $FASTQC "$r"_cutadapt_sickle.fastq.gz \
                --outdir "$curr" \
                -t $NTHREADS &> "$curr"/"$r"_fastqc.log
    done    

    source $VIRTENVS/bwa-meth/bin/activate
    
    fw="$r1"_cutadapt_sickle.fastq.gz
    rv="$r2"_cutadapt_sickle.fastq.gz

    ## bwameth doesn't like reads with different read names, treats them as single end
    # zcat $fw | sed 's/\.1 /\ 1 /g' | gzip -c  > "$fw"_ed.gz
    # zcat $rv | sed 's/\.2 /\ 2 /g'  | gzip -c > "$rv"_ed.gz
    # mv -f "$fw"_ed.gz  "$fw"
    # mv -f "$rv"_ed.gz "$rv"

    bam="$sample"_bwameth_default.bam

    ( bwameth.py --reference "$MM9" \
                 "$fw" "$rv" \
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
             -o mapq.bam
    
    mv -f mapq.bam "$bam"

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
         OUTPUT=$WD/"$(basename $bam .bam)""_dup_removed.bam" \
         METRICS_FILE=$WD/"$(basename $bam .bam)""_dup_removal.metrics"

    mv -f "$bam" "$(basename $bam .bam)"with_duplicates.bam
    mv -f $WD/"$(basename $bam .bam)"_dup_removed.bam "$bam"
    
    # this call does not take into account the min depth nor the mapq thres! methyldackel
    # does not (might not?) filter for that without the mergecontext stuff
    $METHYLDACKEL extract \
                  -q $MAPQ_THRES \
                  -@ $NTHREADS \
                  -d $MIN_DEPTH \
                  --cytosine_report \
                  --CHH \
                  --CHG \
                  -o "$sample"_all_c_"$MIN_DEPTH"_mapq_"$MAPQ_THRES" \
                  $MM9 \
                  "$bam"
    
done < quickchanges.conf


echo cph protocol

while IFS='' read -r line || [[ -n "$line" ]]
do
    cd "$WD"
    
    sample=$(echo $line | cut -f1 -d",")
    r1=$(echo $line | cut -f2 -d",")
    r2=$(echo $line | cut -f3 -d",")

    echo "$(date) cph calls for sample $sample"

done < quickchanges.conf
