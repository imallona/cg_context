#!/bin/bash
##
## Playground with sampled BAM files
##
# Izaskun Mallona
# Mon  9 Apr 08:41:47 CEST 2018

TASK="cg_context"
WD=/home/imallona/"$TASK"

USER="Izaskun Mallona"
HOMEMOUNT=/mnt/home
BAUBECMOUNT=/mnt/baubec
BAUBEC="Group Baubec"

mkdir -p $WD
cd $WD

sudo mount.cifs  //130.60.120.7/home $HOMEMOUNT \
  -o user='Izaskun Mallona',sec=ntlm,uid=1000,gid=1000,iocharset=utf8,file_mode=0777,dir_mode=0777

sudo mount.cifs  //130.60.120.7/"$BAUBEC" $BAUBECMOUNT \
  -o user='Izaskun Mallona',sec=ntlm,uid=1000,gid=1000,iocharset=utf8,file_mode=0777,dir_mode=0777

sudo mount.cifs  //130.60.120.7/"$BAUBEC" /home/imallona/tmp/test \
     -o user='Izaskun Mallona',sec=ntlm,uid=1000,gid=1000,iocharset=utf8,file_mode=0444,dir_mode=044

sudo chmod a+r /home/imallona/tmp/test

## sampling some bamfiles to play with
## TOYBAM=/home/ubuntu/MOUNT/manzo/WGBS/Baubec2015/SRR1653152_trimmed_bismark_bt2.deduplicated.bam

scp -i ~/.ssh/cloudServer.key ubuntu@172.23.28.228:/"$TOYBAM" .

toybam=$(basename "$TOYBAM" .bam)

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

/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump -I --split-files SRR3106762
