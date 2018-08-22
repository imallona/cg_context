#!/bin/bash
##
## methyldackel calls for non cpg methylation
##
# Izaskun Mallona
# aug the 22nd
# GPL

export HOME=/home/imallona
export TASK="cg_context"
export WD="$HOME"/mnt/nfs/"$TASK"
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=16

export MAPQ_THRES=40
export MIN_DEPTH=5

export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel

cd $WD

## samples to call the non cpg methylation from
## here rather than a report we require a minimum depth etc from the very beginning

cat << EOF  >> for_cph.conf
20151223.B-MmES_TKOD3A1c1-3_R,tko+d3a1,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_001,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_001,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_002,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_003,tko+d3a2,bwa_hiseq2k
SRR1274743,tko+d3a2,bwa_miseq_pe
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_001,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_001,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_002,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_003,tko+d3b1,bwa_hiseq2k
SRR1274744,tko+d3b1,bwa_miseq_pe
SRR1274745,tko+d3b1,bwa_miseq_pe
SRR1653150,qko+d3b1,bwa_hiseq2k
SRR1653151,qko+d3b1,bwa_hiseq2k
SRR1653152,qko+d3b1,bwa_hiseq2k
SRR1653153,qko+d3b1,bwa_hiseq2k
SRR1653154,qko+d3b1,bwa_hiseq2k
SRR1653155,qko+d3b1,bwa_hiseq2k
SRR1653156,qko+d3b1,bwa_hiseq2k
SRR1653157,qko+d3b1,bwa_hiseq2k
SRR1653158,qko+d3b1,bwa_hiseq2k
SRR1653159,qko+d3b1,bwa_hiseq2k
SRR1653160,qko+d3b1,bwa_hiseq2k
SRR1653161,qko+d3b1,bwa_hiseq2k
SRR1653162,qko+d3b1,bwa_miseq_pe 
EOF

while IFS='' read -r line || [[ -n "$line" ]]
do
    bam=""
    cd "$WD"

    sample=$(echo $line | cut -f1 -d ',')

    bam=$(find $WD -name "*""$sample""*bwameth_default.bam")
    
    echo $sample
    
    if [ -f $bam ]
    then
        echo 'methyldackel call here'
    else
        echo 'bam not found for ' "$sample"
    fi
    
done < for_cph.conf



# bware the mindepth requires merging the bamfiles first, or representation will get biased
# @todo
# $METHYLDACKEL extract \
#               -q $MAPQ_THRES \
#               -@ $NTHREADS \
#               -d $MIN_DEPTH \
#               --keepStrand \
#               --CHH \
#               -o test.methyldackel \
#               $MM9 \
#               $bam


# awk '
# { 
#   OFS=FS="\t";
#    print $1,$2-1,$2,$4,$5,$3,$7;
# }
# ' "$(basename $bam .bam)".cytosine_report.txt  |
#     "$BEDTOOLS" slop -i - \
#                 -g "$WD"/mm9.genome \
#                 -l 3 -r 4 -s | \
#     "$BEDTOOLS" getfasta -fi $MM9 \
#                 -bed - \
#                 -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
#                 -tab \
#                 -s

# paste "$(basename $bam .bam)".cytosine_report.txt \
#       "$(basename $bam .bam)"_cytosine_report_slop.fa > tmp

# rm -rf "$(basename $bam .bam)"_cytosine_report_slop.fa
# ## now get odd and even lines
# awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tmp > bar
# rm -f tmp 
# mv -f bar "$(basename $bam .bam)"_stranded.txt

# echo "$(date) Processing sample $sample ended"

# gzip "$(basename $bam .bam)"_stranded.txt

# cd "$WD"
# done    

