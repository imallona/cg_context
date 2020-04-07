#!/bin/bash

# tuncay asks for stadler's strand specific readouts with at least 5 reads each strand
# depends uponn stadler_single_end.sh (cg_context and lsh projects for merging)

# Thu Apr 25 14:26:12 CEST 2019


# cat << EOF >> stadler_es.conf
# GSM748786;MouseES_BisSeq_HiSeq;single;SRR299053,SRR299054,SRR299055	
# GSM748787;MouseES_BisSeq_GAIIx;single;SRR299056,SRR299057,SRR299058,SRR299059,SRR299060,SRR299061,SRR299062	
# EOF


export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=16

export MAPQ_THRES=40

export FASTQC=/usr/local/software/FastQC/fastqc
export SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
export FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump


WD="/home/imallona/mnt/nfs/lsh"
BAM="stadler_merged_mapq_40.bam"
# METH=""


cd $WD
$METHYLDACKEL extract \
	      -q $MAPQ_THRES \
	      -@ $NTHREADS \
	      --cytosine_report \
	      $MM9 \
	      $BAM \
	      -o tuncay_$(basename $BAM .bam)
# merge

target=tuncay_stadler_merged_mapq_40.cytosine_report.txt

# awk 'NR%2{printf "%s ",$0;next;}1' $target | \
#     awk '{OFS=FS="\t"; if ($4+$5 >+ 5 && $11+$12 >= 5) print $0 }' > tuncay_cpg_stadler_5_each_strand.txt

awk 'NR%2==0 {print p","$0;} NR%2 {p=$0;}' "$target" | \
    awk '{OFS=FS="\t"; if ($4+$5 >+ 5 && $10+$11 >= 5) print $0}' > tuncay_cpg_stadler_5_each_strand.txt


gzip tuncay_cpg_stadler_5_each_strand.txt

rm -f $target

## same for non cpg

$METHYLDACKEL extract \
	      -q $MAPQ_THRES \
	      -@ $NTHREADS \
	      --cytosine_report \
              --CHG \
              --CHH \
	      $MM9 \
	      $BAM \
	      -o tuncay_no_cpg_$(basename $BAM .bam)


# filter for 5 reads at least

target=tuncay_no_cpg_stadler_merged_mapq_40.cytosine_report.txt

awk '{OFS=FS="\t"; if ($4 > 0 && $4+$5 >= 5 && $6 != "CG") print $0 }' $target > \
    tuncay_non_cpg_stadler_5_cov_and_detected_methylation.txt

gzip tuncay_non_cpg_stadler_5_cov_and_detected_methylation.txt

rm -f $target
