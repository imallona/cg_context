#!/bin/bash
##
## WGBS analysis of SRA se samples
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

cat << EOF >> qko_3b1_se.conf
GSM1545748;qko-dnmt3b1;single;SRR1653150,SRR1653151,SRR1653152,SRR1653153,SRR1653154,SRR1653155,SRR1653156,SRR1653157,SRR1653158,SRR1653159,SRR1653160,SRR1653161
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

	source $VIRTENVS/bwa-meth/bin/activate
	
	( bwameth.py --reference "$MM9" \
                     "$fw" \
                     --threads $NTHREADS | \
                samtools view --threads $NTHREADS -bS - > \
                         "$bam" ) \
            3>&1 1>&2 2>&3 | tee "$sample"_bwameth_default.log

	deactivate

        samtools sort --threads $NTHREADS $bam > sorted
        mv -f sorted $bam
        
	java -jar -XX:ParallelGCThreads=$NTHREADS \
             "$PICARD" MarkDuplicates INPUT="$bam" \
             REMOVE_DUPLICATES=TRUE \
             REMOVE_SEQUENCING_DUPLICATES=TRUE \
             OUTPUT="$(basename $bam .bam)""_dup_removed.bam" \
             METRICS_FILE="$(basename $bam .bam)""_dup_marked.metrics"

        mv -f $bam "$(basename $bam .bam)""_dup_included.bam"
        mv -f "$(basename $bam .bam)""_dup_removed.bam" $bam

        samtools index -b -@ $NTHREADS $bam
        
	"$QUALIMAP"  bamqc \
		     -bam "$bam" \
		     -gd mm9 \
		     -outdir "$(basename $bam .bam)"_qualimap \
		     --java-mem-size=10G \
		     -nt $NTHREADS

	$METHYLDACKEL extract \
                      -q $MAPQ_THRES \
                      -@ $NTHREADS \
                      --cytosine_report \
                      $MM9 \
                      $bam \
                      -o $(basename $bam .bam)
	
	awk '
{ 
  OFS=FS="\t";
   print $1,$2-1,$2,$4,$5,$3,$7;
}
' "$(basename $bam .bam)".cytosine_report.txt  |
            "$BEDTOOLS" slop -i - \
                        -g "$WD"/mm9.genome \
                        -l 3 -r 4 -s | \
            "$BEDTOOLS" getfasta -fi $MM9 \
                        -bed - \
                        -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
                        -tab \
                        -s
        
	paste "$(basename $bam .bam)".cytosine_report.txt \
              "$(basename $bam .bam)"_cytosine_report_slop.fa > tmp

        rm -rf "$(basename $bam .bam)"_cytosine_report_slop.fa
	## now get odd and even lines
	awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tmp > bar
	rm -f tmp 
	mv -f bar "$(basename $bam .bam)"_stranded.txt

	echo "$(date) Processing sample $sample ended"

	gzip "$(basename $bam .bam)"_stranded.txt
        
	cd "$WD"
    done    
done < qko_3b1_se.conf
