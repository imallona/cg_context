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
NTHREADS=20

QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar

cd $WD

## even though this is (probably) not well mapped since 50 nt long and run with bwa-meth bwa-mem default
toy=SRR1653150_1_cutadapt_sickle_bwameth_default.bam
samtools sort $toy --threads $NTHREADS > sorted.bam
mv sorted.bam $toy

"$QUALIMAP"  bamqc \
             -bam "$toy" \
             -gd mm9 \
             -outdir "$(basename $toy .bam)"_qualimap \
             --java-mem-size=10G \
             -nt $NTHREADS
            

## mark duplicates

## beware of the ram
java -jar -XX:ParallelGCThreads=$NTHREADS \
     $MARKDUPLICATES INPUT=$WD/"$toy" \
     OUTPUT=$WD/"$(basename $toy .bam)""_dup_marked.bam" \
     METRICS_FILE=$WD/"$(basename $toy .bam)""_dup_marked.metrics"

## same for the others

for toy in 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default.bam SRR2878513__bwameth_default.bam
do
    samtools sort $toy --threads $NTHREADS > sorted.bam
    mv sorted.bam $toy
    
    "$QUALIMAP"  bamqc \
             -bam "$toy" \
             -gd mm9 \
             -outdir "$(basename $toy .bam)"_qualimap \
             --java-mem-size=40G \
             -nt $NTHREADS


    java -jar -XX:ParallelGCThreads=$NTHREADS \
         $MARKDUPLICATES INPUT=$WD/"$toy" \
         OUTPUT=$WD/"$(basename $toy .bam)""_dup_marked.bam" \
         METRICS_FILE=$WD/"$(basename $toy .bam)""_dup_marked.metrics"
    
done


# split for/rev strands

## this might not be easy for paired end stuff
# https://broadinstitute.github.io/picard/explaian-flags.html
# http://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/

## methyldackel here


