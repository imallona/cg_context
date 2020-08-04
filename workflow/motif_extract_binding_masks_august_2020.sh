#!/bin/bash
##
##
## 12 nov 2019
## Izaskun mallona


## second motif extract without filtering, and with the extract_motifs_frequency_from_bam_binary.sh
##   updated to contain the mincoverage (as well as speed ups)

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



## these come from recap_by_motif_frequencies.sh
cat << EOF  > "$CONFIG_FILE"
20151223.B-MmES_TKOD3A1c1-3_R,tko+d3a1,bwa_hiseq2k
qko+d3b1_merged,qko+d3b1,merged
tko+d3a1_merged,tko+d3a1,merged
tko+d3a2_merged,tko+d3a2,merged
tko+d3b1_merged,tko+d3b1,merged
EOF

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
    # bam='"'"$sample""*bam"'"'
    # find ~/mnt/nfs/cg_context -name $sample -name "*bam"
    bam=$(find -L $WD -type f -name "*bam" | fgrep $sample| head -1)
    echo $bam
done < "$CONFIG_FILE" | wc -l


while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
  
    bam=$(find -L $WD -type f -name "*bam" | fgrep "$sample"| head -1)

    echo $bam
    
    mkdir "$sample"
    cd "$_"
    
    ln -s $bam

    for mask in $(find ../data -name "*bed$")
    do
        ## untested
        bash $EXTRACT_MOTIFS_FREQUENCY_FROM_BAM_BINARY \
             -b $(basename $bam) \
             -t $NTHREADS \
             -c 1 \
             --bedtools $BEDTOOLS \
             --methyldackel $METHYLDACKEL\
             --intervals "$mask"| tee -a with_masks.log
    done            
    
    cd ..
    
done < "$CONFIG_FILE"

