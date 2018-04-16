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

NTHREADS_BOWTIE=2
NTHREADS_BISMARK=2

# baubeclab soft
BISMARK=/home/ubuntu/Soft/Bismark-0.19.0/bismark

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

# echo paired end


for sample in SRR2878513_ 20151223.B-MmES_TKOD3A1c1-3_R
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
done

