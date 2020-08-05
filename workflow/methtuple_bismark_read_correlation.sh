#!/bin/bash
##
## Evaluates the intra-read CpG - CpH correlation/dependence for one sample (paired end)
##
## Izaskun Mallona
## 5 Aug 2020

##
# for instance SRR1274743 TKO_DNMT3A2_bisulfite_miseq paired end

sample=SRR1274743


export HOME=/home/Shared_s3it/imallona/cg_context/nar_review
export TASK=autocorrelation
export WD="$HOME"/"$TASK"

export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=18

export MAPQ_THRES=40

export FASTQC=/usr/local/software/FastQC/fastqc
export SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
export ##CUTADAPT=/usr/local/bin/cutadapt
export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
export MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
export FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump # @todo update

export ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
export ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"

mkdir -p $WD

cd $_


echo 'uncomment to download'
$FASTQDUMP -I --gzip --split-files $sample
export sample="$sample"
echo $sample

# for r in 1 2
# do
#     curr=${WD}/${sample}
#     mkdir -p $curr
    
#     $FASTQC ${WD}/"$sample"_"$r".fastq.gz --outdir ${curr} \
#             -t $NTHREADS &> ${curr}/${sample}_fastqc.log

# done

# source $VIRTENVS/cutadapt/bin/activate

# cutadapt \
#     -j $NTHREADS \
#     -b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
#     -B $ILLUMINA_UNIVERSAL -B $ILLUMINA \
#     -o "$WD"/"${sample}"_1_cutadapt.fastq.gz \
#     -p "$WD"/"${sample}"_2_cutadapt.fastq.gz \
#     "$DATA"/"$sample"_1.fastq.gz "$DATA"/"$sample"_2.fastq.gz &> "$WD"/"$sample"_cutadapt.log

# deactivate


# "$SICKLE" pe \
#           -f "$WD"/"$sample"_1_cutadapt.fastq.gz \
#           -r "$WD"/"$sample"_2_cutadapt.fastq.gz \
#           -o "$WD"/"$sample"_1_cutadapt_sickle.fastq.gz \
#           -p "$WD"/"$sample"_2_cutadapt_sickle.fastq.gz \
#           -t sanger \
#           -s "$WD"/"$sample"_cutadapt_sickle_singles.fastq.gz \
#           -g &> "$WD"/"$sample"_cutadapt_sickle.log


# rm -f "$WD"/"${sample}"_1_cutadapt.fastq.gz "$WD"/"${sample}"_2_cutadapt.fastq.gz

# for r in 1 2
# do
#     curr="$sample"_"$r"_cutadapt_sickle
#     mkdir -p "$WD"/"$curr"
#     $FASTQC "$WD"/${sample}_"$r"_cutadapt_sickle.fastq.gz \
#             --outdir "$WD"/"$curr" \
#             -t $NTHREADS &> ${curr}/${sample}_"$r"_fastqc.log
# done    

# # bismark here


# fw="$sample"_1_cutadapt_sickle.fastq.gz
# rv="$sample"_2_cutadapt_sickle.fastq.gz

# # ## bwameth doesn't like reads with different read names, treats them as single end
# # zcat $fw | sed 's/\.1 /\ 1 /g' | gzip -c  > "$fw"_ed.gz
# # zcat $rv | sed 's/\.2 /\ 2 /g'  | gzip -c > "$rv"_ed.gz
# # mv "$fw"_ed.gz  "$fw"
# # mv "$rv"_ed.gz "$rv"

# ( "$BISMARK" --bowtie1 \
#              --path_to_bowtie $BOWTIE_PATH \
#              --gzip \
#              --parallel $NTHREADS_BISMARK \
#              --genome $(dirname "$MM9") \
#              --1 "$fw" \
#              --2 "$rv") 2>&1 | tee "$WD"/"$sample"_bismark.log


# ## run methtuple to 4 tuples at most, including cpg and cph stuff

