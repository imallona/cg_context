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
NTHREADS=10

QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
BEDTOOLS="$SOFT"/bedtools/bin/bedtools

MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa

MAPQ_THRES=40


cd $WD

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome


echo 'paired end stuff'


for bam in 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam
do
    echo $bam

    echo 'non cpg skipped'
    # --CHG \
    # --CHH \
    $METHYLDACKEL extract \
                  -q $MAPQ_THRES \
                  -@ $NTHREADS \
                  --cytosine_report \
                  $MM9 \
                  $bam \
                  -o $(basename $bam .bam)

    # now fetching the sequence for each c context

    # awk '{OFS=FS="\t"; print $1,$2,$2+1,$4,$5,$3,$7}' 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.cytosine_report.txt  | head > tmp
    
    # "$BEDTOOLS" slop -i tmp \
    #             -g mm9.genome \
    #             -l 2 -r 2 -s | \
    #     "$BEDTOOLS" getfasta -fi $MM9 \
    #                 -bed - \
    #                 -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
    #                 -tab \
    #                 -s

    # paste tmp "$(basename $bam .bam)"_cytosine_report_slop.fa

    awk '
{ OFS=FS="\t";
  if ($3=="+") print $1,$2,$2+1,$4,$5,$3,$7 ;
  else  print $1,$2-2,$2-1,$4,$5,$3,$7 ;

}
' 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.cytosine_report.txt  | head > tmp
        
    "$BEDTOOLS" slop -i tmp \
                -g mm9.genome \
             -l 1 -r 1 | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                 -bed - \
                 -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
                 -tab \
                 -s

    paste tmp "$(basename $bam .bam)"_cytosine_report_slop.fa
    
done


