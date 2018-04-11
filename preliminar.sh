#!/bin/bash
##
## Playground with sampled BAM files
##
# Izaskun Mallona
# Mon  9 Apr 08:41:47 CEST 2018

TASK="cg_context"
WD=/home/imallona/"$TASK"
DATA="$HOME"/data
SOFT="$HOME"/soft

NTHREADS=6

USER="Izaskun Mallona"
BAUBECMOUNT=/home/imallona/mnt/baubec
BAUBEC="Group Baubec"

# taupo soft
FASTQC=/usr/local/software/FastQC/fastqc
SICKLE="$HOME"/sickle/sickle-1.33/sickle

mkdir -p $WD $BAUBECMOUNT
cd $WD

sudo mount.cifs  //130.60.120.7/"$BAUBEC" $BAUBECMOUNT \
  -o user='Izaskun Mallona',sec=ntlm,uid=1000,gid=1000,iocharset=utf8,file_mode=0777,dir_mode=0777


## sampling some bamfiles to play with
## TOYBAM=/home/ubuntu/MOUNT/manzo/WGBS/Baubec2015/SRR1653152_trimmed_bismark_bt2.deduplicated.bam

# scp -i ~/.ssh/cloudServer.key ubuntu@172.23.28.228:/"$TOYBAM" .

# toybam=$(basename "$TOYBAM" .bam)

# MethylDackel mbias \
#              refgenome \
#              $toybam \
#              "$toybam"_test_mbias

# or, rather with some oocytes stuff from GEO
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86297

# let's play with 
# https://www.ncbi.nlm.nih.gov/sra?term=SRX1535160
# for instance, with sra-dump
# SRR3106762	

# /usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump -I --split-files SRR3106762

## using, rather, the real data

cd "$HOME"/data
/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump -I --gzip --split-files SRR1653150
## wtf these are 50 nt long!

# plus those from the synology nas
# "$BAUBECMOUNT"/REPOSITORY/HTS\ data/FGCZ_backup/hiSeq/p2046/HiSeq_20151223_RUN242_DataDelivery


## this oocyte stuff my help too
# https://www.ncbi.nlm.nih.gov/sra?term=SRX1404264
# SRR2878513

/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump -I --gzip --split-files SRR2878513


## fastqc checks

samples="SRR1653150_1 SRR2878513_1 SRR2878513_2 20151223.B-MmES_TKOD3A1c1-3_R1  20151223.B-MmES_TKOD3A1c1-3_R2"
# samples="20151223.B-MmES_TKOD3A1c1-3_R1 20151223.B-MmES_TKOD3A1c1-3_R2"
for sample in ${samples[@]}
do

    curr=${WD}/${sample}
    mkdir -p $curr
    
    $FASTQC ${DATA}/${sample}.fastq.gz --outdir ${curr} \
            -t $NTHREADS &> ${curr}/${sample}_fastqc.log
done

## needed some trimming and adaptor cutting

# no adaptors found for the single end stuff, let's just sickle it
# file:///home/imallona/cg_context/SRR1653150_1/SRR1653150_1_fastqc.html


"$SICKLE" se \
          -c "$DATA"/SRR1653150_1.fastq.gz \
          -g \
          -t sanger \
          -m "$WD"/SRR1653150_1_sickle.fastq.gz &> "${WD}"/SRR1653150_1_sickle.log


## paired end stuff


# currd=$(basename $1)
sample=$(basename $1 .fastq)
echo "Temptative adaptor trimming, just the illumina adaptor"

cutadapt \
    -b $ILLUMINA \
    -o "$TMP"/"${sample}_cutadapted.fastq" \
    -m $MIN_LENGTH \
    --discard-trimmed \
    $(dirname $1)/"$sample"".fastq" &> "${TMP}/${sample}_cutadapt.log" 

echo "Sickle"

"$SICKLE" pe \
          -c "${TMP}/${sample}_cutadapted.fastq" \
          -g \
          -t sanger \
          -m "${TMP}/${sample}_sickle_combined" \
          -s "$TMP"/"$sample"_sickle_singles_trimmed_file &> "${TMP}/${sample}_sickle.log" 

echo "Fastqc"

"$FASTQC" "${TMP}/${sample}_sickle_combined" -outdir "${TMP}" \
                               -t $NTHREADS &> "${TMP}/${sample}_fastq.log"
