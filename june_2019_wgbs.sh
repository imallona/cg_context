#!/bin/bash
##
## new datasets neurodiff
## methyldackel calls for non cpg methylation
##
## Izaskun Mallona
## 14 june 2019
## GPL

export HOME=/home/imallona
export TASK="cg_context"
export WD=/home/Shared_s3it/imallona/"$TASK"/neuro
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=18
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

export PROJECT=p2046
export DATA_DELIVERY=HiSeq2500_20190524_RUN518_o5662_DataDelivery
export CONFIG_FILE=june_2019.conf

mkdir -p $WD

cd $WD

echo Data retrieval

## from fabric
wget --user imallona -e robots=off --ask-password -r --no-parent -nH \
     --cut-dirs=1 \
     -A '20190524.A-TBWGBS*.fastq.gz' \
     --reject='*index.html*' \
     https://fgcz-gstore.uzh.ch/projects/"$PROJECT"/"$DATA_DELIVERY"


cd $WD

echo "Config file"

# IL 12 and IL 18 are two biol replicates from WT DNMT3A2 in TKO
# IL 14 and 19  are two biol replicates from QC DNMT3A2 in TKO
# all samples are different clones, library prepared in paralel and sequenced in parallel
# the only difference is the DNA extraction for 12 and 14 was done in August, while 18 and 19 was done in November

cat << EOF  > "$CONFIG_FILE"
20190524.A-TBWGBS_A_R1.fastq.gz,bwa_hiseq2500_se
20190524.A-TBWGBS_B_R1.fastq.gz,bwa_hiseq2500_se
20190524.A-TBWGBS_K_R1.fastq.gz,bwa_hiseq2500_se
20190524.A-TBWGBS_W_R1.fastq.gz,bwa_hiseq2500_se
EOF

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
              -s "$sample"_cutadapt_sickle_singles.fastq.gz \
              -g &> "$sample"_cutadapt_sickle.log


    rm -f "$r1"_cutadapt.fastq.gz 
    
    curr="$r"_cutadapt_sickle
    mkdir -p "$curr"
    $FASTQC "$r"_cutadapt_sickle.fastq.gz \
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


EOF <<EOF
echo 'cph motifs counts'


mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome


while IFS='' read -r line || [[ -n "$line" ]]
do
    cd $WD
     
    sample=$(echo $line | cut -f1 -d",")
    cd $sample
    
    bam= bam="$sample"_bwameth_default.bam
    ll -h $bam

    
    echo "$(date) cph calls for sample $sample"
    
    # ## mapq over 40
    # was done before
    
    samtools index -@ $NTHREADS "$bam" "$bam".bai

    ## mind the mapq40 filtering was done before
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

    # count instances by motif
    for item in  cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"
    do
        cut -f9 "$item" | sort | uniq -c | sed 's/^ *//' > \
                                               "$sample"_motif_counts_"$item".txt
    done

    rm -f cg_"$sample" ch_"$sample" cg_meth_"$sample" cg_unmeth_"$sample" \
       ch_meth_"$sample" ch_unmeth_"$sample"

    mv tmp_"$sample" "$sample"_raw_report.txt
    rm -f "$sample"_raw_report.txt

    cd $WD
done < "$CONFIG_FILE"


## tests on further stratifying the contexts
cd $WD

for bam in $(find $WD -name "*default.bam")
do
    echo $bam

    bash extract_motifs_frequency_from_bam.sh \
         -b $bam \
         -t $NTHREADS \
         --bedtools $BEDTOOLS \
         --methyldackel $METHYLDACKEL    
done

EOF
