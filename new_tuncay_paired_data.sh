#!/bin/bash
##
## WGBS analysis of SRA pe samples
##
# Izaskun Mallona
# may 3rd 2018
# GPL

export HOME=/home/imallona
export TASK="cg_context_new_tuncay"

export HOME=/home/imallona
export TASK="cg_context_bulk"
export WD="$HOME"/"$TASK"
export DATA="$HOME"/"$TASK"
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=18

export MAPQ_THRES=40

export FASTQC=/usr/local/software/FastQC/fastqc
export SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
export ##CUTADAPT=/usr/local/bin/cutadapt
export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
export FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump 

export ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
export ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"

mkdir -p $WD
cd $_

cat << EOF >> TKODNMT3A2.conf
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_001,BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_R1_001,BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_R2_001
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_001,BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_R1_001,BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_R2_001
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_002,BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_R1_002,BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_R2_002
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_003,BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_R1_003,BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_R2_003
EOF

cat << EOF >> TKODNMT3B1.conf
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_001,BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_R1_001,BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_R2_001
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_001,BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_R1_001,BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_R2_001
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_002,BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_R1_002,BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_R2_002
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_003,BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_R1_003,BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_R2_003
EOF


cd $_

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome

for fn in TKODNMT3A2.conf TKODNMT3B1.conf
do
    while IFS='' read -r line || [[ -n "$line" ]]; do
        echo "Processing: $line"
        sample=$(echo $line | cut -f1 -d",")
        r1=$(echo $line | cut -f2 -d",")
        r2=$(echo $line | cut -f3 -d",")

        ## scp to imlstaupo

        mkdir "$sample"
        cd "$sample"
        
        scp -i ~/.ssh/cloudServer.key ubuntu@172.23.28.228:/~/MOUNT2/tuncay/TEMP/WGBS_FMI/"$r1".fastq.gz .
        scp -i ~/.ssh/cloudServer.key ubuntu@172.23.28.228:/~/MOUNT2/tuncay/TEMP/WGBS_FMI/"$r2".fastq.gz .

        ## process

        for r in "$r1" "$r2"
        do            
            mkdir -p $curr
            
            $FASTQC "$WD"/"$r".fastq.gz --outdir "$curr" \
                    -t $NTHREADS &> ${curr}/"$r"_fastqc.log
        done

        source $VIRTENVS/cutadapt/bin/activate
        
        cutadapt \
            -j $NTHREADS \
            -b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
            -B $ILLUMINA_UNIVERSAL -B $ILLUMINA \
            -o "${sample}"_1_cutadapt.fastq.gz \
            -p "${sample}"_2_cutadapt.fastq.gz \
            "$r1".fastq.gz "$r2".fastq.gz &> "$sample"/"$sample"_cutadapt.log

        deactivate
        
        "$SICKLE" pe \
                  -f "$WD"/"$sample"_1_cutadapt.fastq.gz \
                  -r "$WD"/"$sample"_2_cutadapt.fastq.gz \
                  -o "$WD"/"$sample"_1_cutadapt_sickle.fastq.gz \
                  -p "$WD"/"$sample"_2_cutadapt_sickle.fastq.gz \
                  -t sanger \
                  -s "$WD"/"$sample"_cutadapt_sickle_singles.fastq.gz \
                  -g &> "$WD"/"$sample"_cutadapt_sickle.log

#         ## till here, untested


#         rm -f "$WD"/"${sample}"_1_cutadapt.fastq.gz "$WD"/"${sample}"_2_cutadapt.fastq.gz
        
#         for r in 1 2
#         do
#             curr="$sample"_"$r"_cutadapt_sickle
#             mkdir -p "$WD"/"$curr"
#             $FASTQC "$WD"/${sample}_"$r"_cutadapt_sickle.fastq.gz \
#                     --outdir "$WD"/"$curr" \
#                     -t $NTHREADS &> ${curr}/${sample}_"$r"_fastqc.log
#         done    

#         source $VIRTENVS/bwa-meth/bin/activate
        
#         fw="$sample"_1_cutadapt_sickle.fastq.gz
#         rv="$sample"_2_cutadapt_sickle.fastq.gz

#         ## bwameth doesn't like reads with different read names, treats them as single end
#         zcat $fw | sed 's/\.1 /\ 1 /g' | gzip -c  > "$fw"_ed.gz
#         zcat $rv | sed 's/\.2 /\ 2 /g'  | gzip -c > "$rv"_ed.gz
#         mv "$fw"_ed.gz  "$fw"
#         mv "$rv"_ed.gz "$rv"

#         ( bwameth.py --reference "$MM9" \
#                      $fw $rv \
#                      --threads $NTHREADS |  \
#                 samtools view -bS - | \
#                 samtools sort  - > \
#                          "$sample"_bwameth_default.bam ) \
#             3>&1 1>&2 2>&3 | tee "$sample"_bwameth_default.log

#         deactivate

#         bam="$sample"_bwameth_default.bam

#         "$QUALIMAP"  bamqc \
#                      -bam "$bam" \
#                      -gd mm9 \
#                      -outdir "$(basename $bam .bam)"_qualimap \
#                      --java-mem-size=10G \
#                      -nt $NTHREADS
        
#         java -jar -XX:ParallelGCThreads=$NTHREADS \
#              $MARKDUPLICATES INPUT=$WD/"$bam" \
#              REMOVE_DUPLICATES=TRUE \
#              REMOVE_SEQUENCING_DUPLICATES=TRUE \
#              OUTPUT=$WD/"$(basename $bam .bam)""_dup_marked.bam" \
#              METRICS_FILE=$WD/"$(basename $bam .bam)""_dup_marked.metrics"


#         $METHYLDACKEL extract \
#                       -q $MAPQ_THRES \
#                       -@ $NTHREADS \
#                       --cytosine_report \
#                       $MM9 \
#                       $bam \
#                       -o $(basename $bam .bam)

#         awk '
# { 
#   OFS=FS="\t";
#   if ($3=="+") print $1,$2-1,$2,$4,$5,$3,$7 ;
#   else  print $1,$2-1,$2,$4,$5,$3,$7 ;
# }
# ' "$(basename $bam .bam)".cytosine_report.txt  |
#             "$BEDTOOLS" slop -i - \
#                         -g mm9.genome \
#                         -l 3 -r 4 -s | \
#             "$BEDTOOLS" getfasta -fi $MM9 \
#                         -bed - \
#                         -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
#                         -tab \
#                         -s
#         paste "$(basename $bam .bam)".cytosine_report.txt \
#               "$(basename $bam .bam)"_cytosine_report_slop.fa > tmp

#         ## now get odd and even lines
#         awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tmp > bar
#         rm -f tmp
#         mv -f bar "$(basename $bam .bam)"_stranded.txt
        
#         rm "$r1".fastq.gz
#         rm "$r2".fastq.gz
        
#     done < "$fn"
done

