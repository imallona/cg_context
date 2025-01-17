#!/bin/bash
##
## annotation and filtering of strand-specific meth calls
##
## Izaskun Mallona
## 24 apr 2018
## GPL

echo 'meant to be run at taupo, after fetching the reports from baubec if necessary'
echo 'this is largely duplicated at bowtie_based.sh, which is the latest edit'

export HOME=/home/imallona
export TASK="cytosine_report"

export WD="$HOME"/"$TASK"
export DATA="$HOME"/"$TASK"
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=10

export MAPQ_THRES=40
export MIN_DEPTH=10

export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

cd $WD

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome


# if generated by bowtie, the zero based stuff goes weird
for fn in $(find . -name "*CpG_report*gz")
do
    echo "$fn"

    echo 'filtering min 10 reads, either meth o unmeth'

   #  # this is stupid, must be paired!!
#     zcat $fn | awk '
# { 
#   OFS=FS="\t";
# if ($4+$5>=10) print $0
# }
# ' > filt
    
#     zcat $fn | awk '
# { 
#   OFS=FS="\t";
#   if ($3=="+") print $1,$2,$2+1,$4,$5,$3,$7 ;
#   else  print $1,$2,$2+1,$4,$5,$3,$7 ;
# }
# ' | "$BEDTOOLS" slop -i - \
#                 -g mm9.genome \
#                 -l 3 -r 4 -s | \
#         "$BEDTOOLS" getfasta -fi $MM9 \
#                     -bed - \
#                     -fo "$(basename $fn)"_slop.fa \
#                     -tab \
#                     -s

#     # crashes due to the fact it goes outside the genome boundaries
    
#     paste filt "$(basename $fn)"_slop.fa > tmp

#     ## now get odd and even lines
#     awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tmp > bar
#     rm -f tmp filt
#     mv -f bar "$(basename $fn)"_stranded.txt


    ## rather, paste the outputs and filter first
#     zcat $fn | awk '{printf "%s%s",$0,(NR%2?FS:RS)}' | \
#         awk '
# { 
#   OFS=FS="\t";
# if (($4+$5>=10) && ($10+$11>=10)) print $0
# }' > tmp


    # this is nice
#     zcat $fn | awk '{OFS=FS="\t";printf "%s%s",$0,(NR%2?FS:RS)}' | \
#         awk '
# { 
#   OFS=FS="\t";
# if (($4+$5>=5) && ($11+$12>=5)) print $0
# }' > tmp

#     zcat $fn | head | awk '{OFS=FS="\t";printf "%s%s",$0,(NR%2?FS:RS)}' | \
#         awk '
# { 
#   OFS=FS="\t";
# if (($4+$5>=5) && ($11+$12>=5)) 
#   print $1,$2,$2+1,$4,$5,$3,$7"\n"$8,$9,$9+1,$11,$12,$10,$14
# }' | "$BEDTOOLS" slop -i - \
#                 -g mm9.genome \
#                 -l 3 -r 4 -s | \
#         "$BEDTOOLS" getfasta -fi $MM9 \
#                     -bed - \
#                     -fo "$(basename $fn)"_slop.fa \
#                     -tab \
    #                     -s > foo


#      zcat $fn | head | awk '{OFS=FS="\t";printf "%s%s",$0,(NR%2?FS:RS)}' | \
#         awk '
# { 
#   OFS=FS="\t";
# if (($4+$5>=5) && ($11+$12>=5)) 
#   print $0
# }'
     
#      zcat $fn | head | awk '{OFS=FS="\t";printf "%s%s",$0,(NR%2?FS:RS)}' | \
#         awk '
# { 
#   OFS=FS="\t";
# if (($4+$5>=5) && ($11+$12>=5)) 
#   print $1,$2,$2+1,$3,$4,$5"\n"$8,$9,$9+1,$10,$11,$12
# }' | "$BEDTOOLS" slop -i - \
#                 -g mm9.genome \
#                 -l 3 -r 4 -s | \
#         "$BEDTOOLS" getfasta -fi $MM9 \
#                     -bed - \
#                     -fo "$(basename $fn)"_slop.fa \
#                     -tab \
#                     -s

     zcat $fn | head | awk '{OFS=FS="\t";printf "%s%s",$0,(NR%2?FS:RS)}' | \
        awk '
{ 
  OFS=FS="\t";
if (($4+$5>=5) && ($11+$12>=5)) 
  print $1,$2,$2+1,$4,$5,$3"\n"$8,$9,$9+1,$11,$12,$10
}' | "$BEDTOOLS" slop -i - \
                -g mm9.genome \
                -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$(basename $fn)"_slop.fa \
                    -tab \
                    -s

     # now merging

     rm -f "$(basename $fn)"_slop.fa

     # paste filt "$(basename $fn)"_slop.fa > tmp

     # ## now get odd and even lines
     # awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tmp > bar
     # rm -f tmp filt
     # mv -f bar "$(basename $fn)"_stranded.txt
     
     ##filtering here

done


## in case there were items from methyldackel
