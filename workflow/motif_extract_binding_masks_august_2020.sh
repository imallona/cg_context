#!/bin/bash
##
##
## 04 August 2020
## Izaskun mallona


## second motif extract without filtering, and with the extract_motifs_frequency_from_bam_binary.sh
##   updated to contain the mincoverage and masks (e.g. regions to focus on)

export HOME=/home/imallona
export TASK="cg_context"

export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=20
export MAPQ_THRES=40
export MIN_DEPTH=0

export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export EXTRACT_MOTIFS_FREQUENCY_FROM_BAM_BINARY=~/src/cg_context/extract_motifs_frequency_from_bam_binary.sh

export WD='/home/Shared_s3it/imallona/cg_context/nar_review'

export CONFIG_FILE=with_masks.conf


mkdir -p $WD

cd $WD

# folder bams should be there
ls -l bams

ln -s $SRC/cg_context/data
ls -l data

## these come from recap_by_motif_frequencies.sh
cat << EOF  > "$CONFIG_FILE"
20151223.B-MmES_TKOD3A1c1-3_R,tko+d3a1,bwa_hiseq2k
tko+d3a1_merged,tko+d3a1,merged
tko+d3a2_merged,tko+d3a2,merged
tko+d3b1_merged,tko+d3b1,merged
qko+d3b1_merged,qko+d3b1,merged
EOF

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
    # bam='"'"$sample""*bam"'"'
    # find ~/mnt/nfs/cg_context -name $sample -name "*bam"
    bam=$(find -L $WD -type f -name "*bam" | fgrep $sample| head -1)
    echo $bam
done < "$CONFIG_FILE" | wc -l

# ok, they are four

## test with a bam and mask only start

# mkdir test
# cd $_

# samtools view -bs 0.005 "../bams/tko+d3b1_merged.bam" > sampled.bam

# bam="sampled.bam"
# intervals="../data/NEW_DNMT3A1_sites.bed"
# nthreads=20
# mincov=1
# export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
# export BEDTOOLS="$SOFT"/bedtools/bin/bedtools


# sample=$(basename $bam .bam)
# samtools index -@ $nthreads "$bam" "$bam".bai

# ## mind the mapq40 filtering was done before (bamfile generation)
# $METHYLDACKEL extract \
#               -@ $nthreads \
#               --cytosine_report \
#               --CHH \
#               --CHG \
#               -o "$(basename $bam)"_ch_d \
#               $MM9 \
#               "$bam"


# report="$bam"_ch_d.cytosine_report.txt

# ## if no bedfile with intervals to focus in, proceed
# ##   else, intersect the report with the interval BEDFILE
# if [[ "$intervals" != "None"  ]] ; then
#     echo 'Applying a mask'        

#     ## only grepping elements with data (not 0,0 for meth and unmeth reads)
#     awk '{FS=OFS="\t"; print $1,$2,$2+1,$3,$4,$5,$6,$7}' "$report" |
#         grep -vP '\t0\t0\t' | \
#         "$BEDTOOLS" intersect -a - -b "$intervals" -wa > masked
    
#     report="$bam"_ch_d.cytosine_report_mask_"$(basename $intervals)".txt
#     mv -f masked "$report"
  
# fi



# bash $EXTRACT_MOTIFS_FREQUENCY_FROM_BAM_BINARY \
#              -b  $bam \
#              -t $NTHREADS \
#              -c 1 \
#              --bedtools $BEDTOOLS \
#              --methyldackel $METHYLDACKEL \
#              --intervals "$mask"| tee -a with_mask_"$(basename $mask .bed)".log



## test with a bam and mask only end


cd $WD

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
  
    bam=$(find -L $WD -type f -name "*bam" | fgrep "$sample"| head -1)

    echo $bam

    mkdir -p "$WD"/"$sample"
    cd "$_"
    
    # ln -s $bam

    for interval in $(find $WD/data/ -name "*bed")
    do
        echo $interval
        
        mkdir -p "$(basename $interval .bed)"
        cd "$_"
        
        ## untested
        bash $EXTRACT_MOTIFS_FREQUENCY_FROM_BAM_BINARY \
             -b $bam \
             -t $NTHREADS \
             -c 1 \
             --bedtools $BEDTOOLS \
             --methyldackel $METHYLDACKEL \
             --intervals "$interval"| tee -a with_mask_"$(basename $interval)".log
    done            
    
    cd $WD
    
done < "$CONFIG_FILE"

