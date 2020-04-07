#!/bin/bash
##
## WGBS analysis of SRA se samples (this is a companion script to download
#  some data while the main tuncay_qko_3b1_se.sh runs the main stuff
##
# Izaskun Mallona
# aug the 14th
# GPL


export HOME=/home/imallona
export TASK="cg_context"
export WD="$HOME"/mnt/nfs/"$TASK"
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=16

export MAPQ_THRES=40

export FASTQC=/usr/local/software/FastQC/fastqc
export SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export PICARD="$SOFT"/picard/build/libs/picard.jar

export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
export FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump 

export ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
export ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"

mkdir -p $WD
cd $_

mysql --user=genome \
      --host=genome-euro-mysql.soe.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > \
      mm9.genome

cat << EOF >> qko_3b1_se_companion.conf
GSM1545748;qko-dnmt3b1;single;SRR1653161,SRR1653160,SRR1653159,SRR1653158,SRR1653157
EOF

while IFS='' read -r line || [[ -n "$line" ]]
do

    cd "$WD"

    samples=$(echo $line | cut -f4 -d ';')
    for sample in  $(echo $samples | tr "," "\n")
    do
	
	echo "$(date) Processing sample $sample"
	## scp to imlstaupo

	mkdir - "$sample"
	cd "$sample"
        
	r="$sample"_1
        fw="$r"_cutadapt_sickle.fastq.gz
        bam="$sample"_bwameth_default.bam
        

        if [ ! -f "$fw" ]
        then
	    $FASTQDUMP -I --gzip --split-files $sample

	    curr="$r"_raw
	    mkdir -p $curr
	    
	    $FASTQC "$r".fastq.gz --outdir "$curr" \
		    -t $NTHREADS &> "$curr"/"$r"_fastqc.log

	    source $VIRTENVS/cutadapt/bin/activate
	    
	    cutadapt \
	        -j $NTHREADS \
	        -b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
	        -o "$r"_cutadapt.fastq.gz \
	        "$r".fastq.gz &> "$sample"_cutadapt.log

	    deactivate

	    rm -f "$r".fastq.gz 
	    
	    "$SICKLE" se \
		      -f "$r"_cutadapt.fastq.gz \
		      -o "$r"_cutadapt_sickle.fastq.gz \
		      -t sanger \
		      -g &> "$sample"_cutadapt_sickle.log


	    rm -f "$r"_cutadapt.fastq.gz
	    
	    curr="$r"_cutadapt_sickle
	    mkdir -p "$curr"
	    $FASTQC "$r"_cutadapt_sickle.fastq.gz \
		    --outdir "$curr" \
		    -t $NTHREADS &> "$curr"/"$r"_fastqc.log
        fi
    done    
done < qko_3b1_se_companion.conf
