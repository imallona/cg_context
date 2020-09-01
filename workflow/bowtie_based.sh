#!/bin/bash
##
## WGBS mapping using bowtie
##
# Izaskun Mallona
# 16th apr 2018
##
## to be run at baubeclab
# GPL

TASK="cg_context"
WD=/home/ubuntu/MOUNT2/imallona/"$TASK"
# DATA="$HOME"/data
DATA=/home/ubuntu/MOUNT2/imallona
SOFT="$HOME"/soft
MM9=/home/ubuntu/REPOS/annotation/mm9_genome/
VIRTENVS=~/virtenvs

NTHREADS_BOWTIE=4
NTHREADS_BISMARK=2
NTHREADS=6

# baubeclab soft
BISMARK=/home/ubuntu/Soft/Bismark-0.19.0/bismark
QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap

MAPQ=40

cd $DATA

echo 'rsync -a -v imallona@imlstaupo.uzh.ch:/home/imallona/cg_context/*cutadapt_sickle.fastq.gz'

mkdir -p $WD
cd $WD

## let's index the genome (rather not, use baubeclab's)
echo single end

for sample in SRR1653150_1
do
    echo $sample
    fastq="$DATA"/"$sample"_cutadapt_sickle.fastq.gz
    ( "$BISMARK" --bowtie2 \
               --genome $MM9 \
               --N 0 \
               -D 20 \
               -R 3 \
               -L 20 \
               --gzip \
               --parallel $NTHREADS_BISMARK \
               -p $NTHREADS_BOWTIE \
               --se "$fastq" ) 2>&1 | tee "$WD"/"$sample"_bismark.log
done

# echo single end


for sample in 20151223.B-MmES_TKOD3A1c1-3_R SRR2878513_
do
    echo $sample
    r1="$DATA"/"$sample"1_cutadapt_sickle.fastq.gz
    r2="$DATA"/"$sample"2_cutadapt_sickle.fastq.gz

    echo $sample
    fastq="$DATA"/"$sample"_cutadapt_sickle.fastq.gz
    ( "$BISMARK" --bowtie2 \
                 --genome $MM9 \
                 --N 0 \
                 -D 20 \
                 -R 3 \
                 -L 20 \
                 --gzip \
                 --parallel $NTHREADS_BISMARK \
                 -p $NTHREADS_BOWTIE \
                 --1 "$r1" \
                 --2 "$r2"
    ) 2>&1 | tee "$WD"/"$sample"_bismark.log

    echo qualimap
    bam="$sample"1_cutadapt_sickle_bismark_bt2_pe.bam
    samtools sort --threads $NTHREADS "$bam" > sorted.bam
    mv -v sorted.bam $bam
    
    "$QUALIMAP"  bamqc \
                 -bam "$bam" \
                 -gd mm9 \
                 -outdir "$(basename $bam .bam)"_qualimap \
                 --java-mem-size=10G \
                 -nt $NTHREADS
    
    echo dedup

    "$(dirname $BISMARK)"/deduplicate_bismark -p \
                         --output_dir $WD \
                         --bam \
                         "$sample"1_cutadapt_sickle_bismark_bt2_pe.bam


    echo meth extract

    # bismark_methylation_extractor --bedGraph --gzip Sample2_SE_trimmed.fq.gz_bismark.bam
    "$(dirname $BISMARK)"/bismark_methylation_extractor \
                         --paired \
                         --parallel $NTHREADS \
                         --zero_based \
                         --cytosine_report \
                         --genome_folder $MM9 \
                         "$sample"1_cutadapt_sickle_bismark_bt2_pe.bam
    
done


## broken at some point, let-s call the meth since de dedup worked
bam=SRR1653150_1_cutadapt_sickle_bismark_bt2.bam
"$(dirname $BISMARK)"/deduplicate_bismark -s \
                     --output_dir $WD \
                     --bam \
                     "$bam"


# for sample in 20151223.B-MmES_TKOD3A1c1-3 SRR2878513
# do
   
#     echo meth extract

#     "$(dirname $BISMARK)"/bismark_methylation_extractor \
#                          --paired \
#                          --parallel $NTHREADS \
#                          --zero_based \
#                          --cytosine_report \
#                          --genome_folder $MM9 \
#                          "$sample"*_cutadapt_sickle_bismark_bt2_*e.deduplicated.bam
    
# done


# ## again, deduplicating 150

# for se_sample in SRR1653150
# do
   
#     echo meth extract
#     echo "$se_sample"

#     "$(dirname $BISMARK)"/bismark_methylation_extractor \
#                          --single-end \
#                          --parallel $NTHREADS \
#                          --zero_based \
#                          --cytosine_report \
#                          --genome_folder $MM9 \
#                          "$se_sample"*_cutadapt_sickle_bismark_bt2.deduplicated.bam
    
# done



for sample in 20151223.B-MmES_TKOD3A1c1-3 SRR2878513
do
   
    echo meth extract

    "$(dirname $BISMARK)"/bismark_methylation_extractor \
                         --paired \
                         --parallel $NTHREADS \
                         --zero_based \
                         --cytosine_report \
                         --genome_folder $MM9 \
                         "$sample"*_cutadapt_sickle_bismark_bt2_*e.deduplicated.bam
    
done


## again, calls after mapq filtering
# cannot pipe this properly

for se_sample in SRR1653150
do
   
    echo meth extract
    echo "$se_sample"

    samtools view \
             -q $MAPQ \
             -@ $NTHREADS \
             -b "$se_sample"*_cutadapt_sickle_bismark_bt2.deduplicated.bam > \
             "$se_sample"*_cutadapt_sickle_bismark_bt2.deduplicated_mapq"$MAPQ".bam

    "$(dirname $BISMARK)"/bismark_methylation_extractor \
                         --single-end \
                         --parallel $NTHREADS \
                         --zero_based \
                         --cytosine_report \
                         --genome_folder $MM9 \
                         --gzip \
                         "$se_sample"*_cutadapt_sickle_bismark_bt2.deduplicated_mapq"$MAPQ".bam

    rm  "$se_sample"*_cutadapt_sickle_bismark_bt2.deduplicated_mapq"$MAPQ".bam
    
done

for sample in 20151223.B-MmES_TKOD3A1c1-3 SRR2878513
do
   
    echo meth extract
    echo $sample

    samtools view \
             -@ $NTHREADS \
             -q $MAPQ \
             -b  "$sample"*_cutadapt_sickle_bismark_bt2_*e.deduplicated.bam > \
             "$sample"*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq"$MAPQ".bam

    "$(dirname $BISMARK)"/bismark_methylation_extractor \
                         --paired \
                         --parallel $NTHREADS \
                         --zero_based \
                         --cytosine_report \
                         --gzip \
                         --genome_folder $MM9 \
                         "$sample"*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq"$MAPQ".bam

    rm "$sample"*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq"$MAPQ".bam
    
done


## strand determinations




## this rather at taupo

export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa

export HOME=/home/imallona
export TASK="cg_context_bulk"

export WD="$HOME"/"$TASK"
export DATA="$HOME"/"$TASK"
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export NTHREADS=18


mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome

ln -s /home/imallona/mnt/baubec/cg_context_from_baubec/20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40.CpG_report.txt.CpG_report.txt.gz


# if generated by bowtie, the zero based stuff goes weird
# for fn in $(find . -name "*CpG_report*gz")
# for fn in "./SRR1653150*_cutadapt_sickle_bismark_bt2.deduplicated_mapq40.CpG_report.txt.CpG_report.txt.gz"
for fn in 20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40.CpG_report.txt.CpG_report.txt.gz
do
    echo "$fn"

    zcat "$fn" | awk '{OFS=FS="\t";printf "%s%s",$0,(NR%2?FS:RS)}' | \
        awk '
    { 
      OFS=FS="\t";
      print $1,$2,$2+1,$4,$5,$3"\n"$8,$9,$9+1,$11,$12,$10
    }' | "$BEDTOOLS" slop -i - \
                     -g mm9.genome \
                     -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$(basename "$fn")"_slop.fa \
                    -tab \
                    -s

    zcat $fn | \
        paste - "$(basename "$fn")"_slop.fa | \
        awk '{printf "%s%s",$0,(NR%2?FS:RS)}' > \
            "$(basename $fn .CpG_report.txt.CpG_report.txt.gz)"_stranded.txt
    
    gzip "$(basename $fn .CpG_report.txt.CpG_report.txt.gz)"_stranded.txt
done
