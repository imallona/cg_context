#!/bin/bash
##
## BAM files parsing to call strand specific methylation
##
# Izaskun Mallona
# Mon  12 apr 2018
# GPL
##
## requires preliminar.sh to be run first

TASK="cg_context"
WD=/home/imallona/"$TASK"
DATA="$HOME"/data
SOFT="$HOME"/soft
NTHREADS=6

QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel

cd $WD

## even though this is (probably) not well mapped since 50 nt long and run with bwa-meth bwa-mem default
toy=SRR1653150_1_cutadapt_sickle_bwameth_default.bam
samtools sort $toy --threads $NTHREADS > sorted.bam
mv sorted.bam SRR1653150_1_cutadapt_sickle_bwameth_default_sorted.bam

"$QUALIMAP"  bamqc \
             -bam SRR1653150_1_cutadapt_sickle_bwameth_default_sorted.bam \
             -gd mm9 \
             -outdir "$(basename $toy .bam)"_qualimap \
             --java-mem-size=10G \
             -nt $NTHREADS
            

## mark duplicates

## methyldackel here
