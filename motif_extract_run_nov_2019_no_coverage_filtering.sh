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

export WD='/home/imallona/cg_context/nov_2019_no_coverage_filtering'

export CONFIG_FILE=nov_2019_no_coverage_filtering.conf


mkdir -p $WD

cd $WD



## these come from recap_by_motif_frequencies.sh
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

## and these from june_2019_wgbs.sh

cat << EOF  >> "$CONFIG_FILE"
20190524.A-TBWGBS_A_R1.fastq.gz,bwa_hiseq2500_se
20190524.A-TBWGBS_B_R1.fastq.gz,bwa_hiseq2500_se
20190524.A-TBWGBS_K_R1.fastq.gz,bwa_hiseq2500_se
20190524.A-TBWGBS_W_R1.fastq.gz,bwa_hiseq2500_se
EOF



# old location
ln -s ~/mnt/nfs/cg_context
ln -s /home/Shared_s3it/imallona/cg_context/neuro

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
    # bam='"'"$sample""*bam"'"'
    # find ~/mnt/nfs/cg_context -name $sample -name "*bam"
    bam=$(find -L $WD -type f -name "*bam" | fgrep $sample| head -1)
    echo $bam
done < "$CONFIG_FILE" | wc -l

# 39 as well, cool



while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
  
    bam=$(find -L $WD -type f -name "*bam" | fgrep "$sample"| head -1)

    echo $bam
    
    mkdir "$sample"
    cd "$_"
    
    ln -s $bam

    bash $EXTRACT_MOTIFS_FREQUENCY_FROM_BAM_BINARY \
         -b $(basename $bam) \
         -t $NTHREADS \
         -c 0 \
         --bedtools $BEDTOOLS \
         --methyldackel $METHYLDACKEL | tee -a no_filtering.log

    cd ..
    
done < "$CONFIG_FILE"

