#!/bin/bash
##
## everything run with proportions
##
# Izaskun Mallona
## 21feb 2019
# GPL

export HOME=/home/imallona
export TASK="cg_context"
export WD="$HOME"/mnt/nfs/"$TASK"/feb_2019
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=12
export MAPQ_THRES=40
export MIN_DEPTH=10

export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export CONFIG_FILE=feb_2019.conf

mkdir -p $WD

cd $WD

echo "Config file"

# merged datasets come from `cph_calls.sh` and mix seq platforms
# 20151223.B-MmES_TKOD3A1c1-3_R,tko+d3a1,bwa_hiseq2k
# BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_001,tko+d3a2,bwa_hiseq2k
# BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_001,tko+d3a2,bwa_hiseq2k
# BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_002,tko+d3a2,bwa_hiseq2k
# BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_003,tko+d3a2,bwa_hiseq2k
# SRR1274743,tko+d3a2,bwa_miseq_pe
# BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_001,tko+d3b1,bwa_hiseq2k
# BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_001,tko+d3b1,bwa_hiseq2k
# BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_002,tko+d3b1,bwa_hiseq2k
# BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_003,tko+d3b1,bwa_hiseq2k
# SRR1274744,tko+d3b1,bwa_miseq_pe
# SRR1274745,tko+d3b1,bwa_miseq_pe
# SRR1653150,qko+d3b1,bwa_hiseq2k
# SRR1653151,qko+d3b1,bwa_hiseq2k
# SRR1653152,qko+d3b1,bwa_hiseq2k
# SRR1653153,qko+d3b1,bwa_hiseq2k
# SRR1653154,qko+d3b1,bwa_hiseq2k
# SRR1653155,qko+d3b1,bwa_hiseq2k
# SRR1653156,qko+d3b1,bwa_hiseq2k
# SRR1653157,qko+d3b1,bwa_hiseq2k
# SRR1653158,qko+d3b1,bwa_hiseq2k
# SRR1653159,qko+d3b1,bwa_hiseq2k
# SRR1653160,qko+d3b1,bwa_hiseq2k
# SRR1653161,qko+d3b1,bwa_hiseq2k
# SRR1653162,qko+d3b1,bwa_miseq_pe 


cat << EOF  > "$CONFIG_FILE"
20181207.A-WGBS_IL12,20181207.A-WGBS_IL12_R1,20181207.A-WGBS_IL12_R2,bwa_hiseq4k_pe
20181207.A-WGBS_IL14,20181207.A-WGBS_IL14_R1,20181207.A-WGBS_IL14_R2,bwa_hiseq4k_pe
20181207.A-WGBS_IL18,20181207.A-WGBS_IL18_R1,20181207.A-WGBS_IL18_R2,bwa_hiseq4k_pe
20181207.A-WGBS_IL19,20181207.A-WGBS_IL19_R1,20181207.A-WGBS_IL19_R2,bwa_hiseq4k_pe
20180913.A-WGBS_3,20180913.A-WGBS_3_R1,20180913.A-WGBS_3_R2,bwa_hiseq4k_pe
20180913.A-WGBS_4,20180913.A-WGBS_4_R1,20180913.A-WGBS_4_R2,bwa_hiseq4k_pe
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
qko+d3b1_merged,qko+d3b1,merged
tko+d3a1_merged,tko+d3a1,merged
tko+d3a2_merged,tko+d3a2,merged
tko+d3b1_merged,tko+d3b1,merged
EOF

## old datasets, merged and without merging



while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
    # bam='"'"$sample""*bam"'"'
    # find ~/mnt/nfs/cg_context -name $sample -name "*bam"
    bam=$(find ~/mnt/nfs/cg_context -name "*bam" | fgrep $sample| head -1)
    echo $bam
done < "$CONFIG_FILE" | wc -l

wc -l $CONFIG_FILE

## ok, the 35 files found

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
  
    bam=$(find ~/mnt/nfs/cg_context -name "*bam" | fgrep "$sample"| head -1)

    mkdir "$sample"
    cd "$_"

    bash $SRC/cg_context/extract_motifs_frequency_from_bam.sh \
         -b "$bam" \
         -t $NTHREADS \
         --bedtools $BEDTOOLS \
         --methyldackel $METHYLDACKEL

    cd ..
    
done < "$CONFIG_FILE"
