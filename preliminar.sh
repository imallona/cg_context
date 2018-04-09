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
