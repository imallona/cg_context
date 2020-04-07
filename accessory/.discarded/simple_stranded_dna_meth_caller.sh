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


# for bam in 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam SRR2878513__bwameth_default_dup_marked.bam
for bam in SRR2878513__bwameth_default_dup_marked.bam
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

#     awk '
# { OFS=FS="\t";
#   if ($3=="+") print $1,$2,$2+1,$4,$5,$3,$7 ;
#   else  print $1,$2-2,$2-1,$4,$5,$3,$7 ;

# }
# ' "$(basename $bam .bam)".cytosine_report.txt  |
#     "$BEDTOOLS" slop -i - \
#                 -g mm9.genome \
#                 -l 3 -r 3 | \
#     "$BEDTOOLS" getfasta -fi $MM9 \
#                 -bed - \
#                 -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
#                 -tab \
#                 -s

#     awk '
# { OFS=FS="\t";
#   if ($3=="+") print $1,$2-4,$2+4,$4,$5,$3,$7 ;
#   else  print $1,$2-5,$2+3,$4,$5,$3,$7 ;

# }
# ' "$(basename $bam .bam)".cytosine_report.txt  |
#         "$BEDTOOLS" slop -i - \
#                     -g mm9.genome \
#                     -l 0 -r 0 | \
#         "$BEDTOOLS" getfasta -fi $MM9 \
#                     -bed - \
#                     -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
#                     -tab \
#                     -s

    #     # crashes due to the fact it goes outside the genome boundaries

    awk '
{ OFS=FS="\t";
  if ($3=="+") print $1,$2-1,$2,$4,$5,$3,$7 ;
  else  print $1,$2-1,$2,$4,$5,$3,$7 ;

}
' "$(basename $bam .bam)".cytosine_report.txt  |
        "$BEDTOOLS" slop -i - \
                    -g mm9.genome \
                    -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
                    -tab \
                    -s

    # crashes due to the fact it goes outside the genome boundaries
    
    paste "$(basename $bam .bam)".cytosine_report.txt \
          "$(basename $bam .bam)"_cytosine_report_slop.fa > tmp

    ## now get odd and even lines
    awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tmp > bar
    rm -f tmp
    mv bar "$(basename $bam .bam)"_stranded.txt
    


done


