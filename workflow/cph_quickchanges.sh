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
export MAPQ_THRES=40

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
    
    # java -jar -XX:ParallelGCThreads=$NTHREADS \
    #      $MARKDUPLICATES INPUT=$WD/"$bam" \
    #      REMOVE_DUPLICATES=TRUE \
    #      REMOVE_SEQUENCING_DUPLICATES=TRUE \
    #      OUTPUT=$WD/"$(basename $bam .bam)""_dup_removed.bam" \
    #      METRICS_FILE=$WD/"$(basename $bam .bam)""_dup_removal.metrics"

    java -jar -XX:ParallelGCThreads=$NTHREADS \
     "$PICARD" MarkDuplicates INPUT=$WD/"$bam" \
     REMOVE_DUPLICATES=TRUE \
     REMOVE_SEQUENCING_DUPLICATES=TRUE \
     OUTPUT=$WD/"$(basename $bam .bam)""_dup_removed.bam" \
     METRICS_FILE=$WD/"$(basename $bam .bam)""_dup_marked.metrics"

    mv -f "$bam" "$(basename $bam .bam)"with_duplicates.bam
    mv -f $WD/"$(basename $bam .bam)"_dup_removed.bam "$bam"
    
    # # this call does not take into account the min depth nor the mapq thres! methyldackel
    # # does not (might not?) filter for that without the mergecontext stuff
    # $METHYLDACKEL extract \
    #               -q $MAPQ_THRES \
    #               -@ $NTHREADS \
    #               -d $MIN_DEPTH \
    #               --cytosine_report \
    #               --CHH \
    #               --CHG \
    #               -o "$sample"_all_c_"$MIN_DEPTH"_mapq_"$MAPQ_THRES" \
    #               $MM9 \
    #               "$bam"
    
done < quickchanges.conf


echo cph protocol

while IFS='' read -r line || [[ -n "$line" ]]
do
    cd "$WD"
    
    sample=$(echo $line | cut -f1 -d",")
    bam="$sample"_bwameth_default.bam

    echo "$(date) cph calls for sample $sample"
    
    ## mapq over 40
    samtools view -@ $NTHREADS -h -b -q "$MAPQ_THRES" "$bam" \
             -o "$sample"_mapq.bam
    
    mv -f "$sample"_mapq.bam "$bam"
    
    samtools index -@ $NTHREADS "$bam" "$bam".bai

    # this call does not take into account the min depth nor the mapq thres! methyldackel
    # does not (might not?) filter for that without the mergecontext stuff
    $METHYLDACKEL extract \
                  -q $MAPQ_THRES \
                  -@ $NTHREADS \
                  -d $MIN_DEPTH \
                  --cytosine_report \
                  --CHH \
                  --CHG \
                  -o "$bam"_ch_d_"$MIN_DEPTH"_mapq_"$MAPQ_THRES" \
                  $MM9 \
                  "$bam"

    report="$bam"_ch_d_"$MIN_DEPTH"_mapq_"$MAPQ_THRES".cytosine_report.txt

    # only stuff  with a coverage of at least 10
    # echo 'getting 10k covered Cs only'
    awk '{
FS=OFS="\t"; 
if ($4+$5 >= 10)
 print $0
}' \
        "$report" > tmp_covered_"$sample"

    mv -f tmp_covered_"$sample" "$report"
    
    ## motif retrieval

    awk '
{ 
  OFS=FS="\t";
   print $1,$2-1,$2,$4,$5,$3,$7;
}
' "$report"  |
        "$BEDTOOLS" slop -i - \
                    -g "$WD"/mm9.genome \
                    -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$report"_cytosine_report_slop.fa \
                    -tab \
                    -s
    paste "$report" \
          "$report"_cytosine_report_slop.fa > tmp_"$sample"

    rm -rf "$report"_cytosine_report_slop.fa
    # get meth and unmeth stats
    
    # split into cg and non cg 
    awk '{OFS=FS="\t"; if ($6 == "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' \
        tmp_"$sample" > cg_"$sample"
    awk '{OFS=FS="\t"; if ($6 != "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' \
        tmp_"$sample" > ch_"$sample"
    
    # split into meth and unmeth
    awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' cg_"$sample" > cg_meth_"$sample"
    awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' cg_"$sample" > cg_unmeth_"$sample"
    
    awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' ch_"$sample" > ch_meth_"$sample"
    awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' ch_"$sample" > ch_unmeth_"$sample"
    wc -l cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"

    # count instancees by motif
    for item in  cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"
    do
        cut -f9 "$item" | sort | uniq -c | sed 's/^ *//' > \
                                               "$sample"_motif_counts_"$item".txt
    done

    rm -f cg_"$sample" ch_"$sample" cg_meth_"$sample" cg_unmeth_"$sample" \
       ch_meth_"$sample" ch_unmeth_"$sample"

    mv tmp_"$sample" "$sample"_raw_report.txt
    gzip "$sample"_raw_report.txt

done < quickchanges.conf
