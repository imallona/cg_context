#!/bin/bash
##
## Playground with sampled BAM files
##
# Izaskun Mallona
# Mon  9 Apr 08:41:47 CEST 2018

TASK="cg_context"
WD=/home/imallona/"$TASK"

USER="Izaskun Mallona"

BAUBECMOUNT=/home/imallona/mnt/baubec
BAUBEC="Group Baubec"

# taupo
FASTQC=

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

# for sample in ${samples[@]}
# do
#     for pair in R1 R2
#     do
#         curr=${WD}/quality/${sample}_${pair}
#         mkdir $curr
        
#         $FASTQC ${DAT}/${sample}_${pair}.fastq --outdir ${curr} \
#                 -t $NTHREADS &> ${curr}/${sample}_${pair}_fastqc.log
#     done
# done
