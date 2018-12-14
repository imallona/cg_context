#!/bin/bash
##
## new datasets
## methyldackel calls for non cpg methylation
##
# Izaskun Mallona
## 14th dec 2018
# GPL

export HOME=/home/imallona
export TASK="cg_context"
export WD="$HOME"/mnt/nfs/"$TASK"/dec_2018
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=18
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

export PROJECT=p2046
export DATA_DELIVERY=HiSeq4000_20181207_RUN499_o5114_DataDelivery
export CONFIG_FILE=dec_2018.conf

mkdir -p $WD

cd $WD

echo Data retrieval

## from fabric
wget --user imallona -e robots=off --ask-password -r --no-parent -nH \
     --cut-dirs=1 \
     -A '*WGBS*.fastq.gz' \
     --reject='*index.html*' \
     https://fgcz-gstore.uzh.ch/projects/"$PROJECT"/"$DATA_DELIVERY"


cd $WD

echo "Config file"

# IL 12 and IL 18 are two biol replicates from WT DNMT3A2 in TKO
# IL 14 and 19  are two biol replicates from QC DNMT3A2 in TKO
# all samples are different clones, library prepared in paralel and sequenced in parallel
# the only difference is the DNA extraction for 12 and 14 was done in August, while 18 and 19 was done in November

cat << EOF  > "$CONFIG_FILE"
20181207.A-WGBS_IL12,20181207.A-WGBS_IL12_R1,20181207.A-WGBS_IL12_R2,bwa_hiseq4k_pe
20181207.A-WGBS_IL14,20181207.A-WGBS_IL14_R1,20181207.A-WGBS_IL14_R2,bwa_hiseq4k_pe
20181207.A-WGBS_IL18,20181207.A-WGBS_IL18_R1,20181207.A-WGBS_IL18_R2,bwa_hiseq4k_pe
20181207.A-WGBS_IL19,20181207.A-WGBS_IL19_R1,20181207.A-WGBS_IL19_R2,bwa_hiseq4k_pe
EOF

echo Mapping


while IFS='' read -r line || [[ -n "$line" ]]; do
    cd "$WD"
    
    sample=$(echo $line | cut -f1 -d",")
    r1=$(echo $line | cut -f2 -d",")
    r2=$(echo $line | cut -f3 -d",")

    echo "$(date) Processing sample $sample"
    echo "... read 1 is $r1"
    echo "... read 2 is $r2"

    ## process
    mkdir -p "$sample"
    cd $_
    
    for r in "$r1" "$r2"
    do
        ln -s "$WD"/"$PROJECT"/"$DATA_DELIVERY"/"$r".fastq.gz
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
             -o "$sample"_mapq.bam
    
    mv -f "$sample"_mapq.bam "$bam"

    "$QUALIMAP"  bamqc \
                 -bam "$bam" \
                 -gd mm9 \
                 -outdir "$(basename $bam .bam)"_qualimap \
                 --java-mem-size=10G \
                 -nt $NTHREADS

    java -jar -XX:ParallelGCThreads=$NTHREADS \
     "$PICARD" MarkDuplicates INPUT="$bam" \
     REMOVE_DUPLICATES=TRUE \
     REMOVE_SEQUENCING_DUPLICATES=TRUE \
     OUTPUT="$(basename $bam .bam)""_dup_removed.bam" \
     METRICS_FILE="$(basename $bam .bam)""_dup_marked.metrics"

    mv -f "$bam" "$(basename $bam .bam)"with_duplicates.bam
    mv -f "$(basename $bam .bam)"_dup_removed.bam "$bam"

    cd $WD
        
done < "$CONFIG_FILE"






EOF <<EOF

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome




## bulk bamfiles check

for fn in $(find .. -name "*bam" ! -name "*dup_included*" ! -name "*with_duplicates*" \
                  ! -name "*merged*")
do
    bam=$(basename $fn)
    sample=$(basename $bam .bam)
    echo $sample    
    ###
    
    echo "$(date) cph calls for sample $sample"
    
    # ## mapq over 40
    # samtools view -@ $NTHREADS -h -b -q "$MAPQ_THRES" "$bam" \
    #          -o "$sample"_mapq.bam
    
    # mv -f "$sample"_mapq.bam "$bam"

    ln -s $fn
    
    samtools index -@ $NTHREADS "$bam" "$bam".bai

    # this call does not take into account the min depth nor the mapq thres! methyldackel
    # does not (might not?) filter for that without the mergecontext flag
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
    rm -f "$sample"_raw_report.txt

done 

## cleaning (removal of .bai s etc)


EOF
