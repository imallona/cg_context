#!/bin/bash
#
# deptools based bam file comparison
# Apr 26 2018
#
# Izaskun Mallona
#
# GPLv2


export HOME=/home/imallona
export TASK="cg_context_bulk"

export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=10

export MAPQ_THRES=40

# export FASTQC=/usr/local/software/FastQC/fastqc
# export SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
# export ##CUTADAPT=/usr/local/bin/cutadapt
# export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
# export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
# export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
# export MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
# export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
# export FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump 

# export ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
# export ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"


source "$VIRTENVS"/deeptools/bin/activate


multiBamSummary bins \
                --numberOfProcessors $NTHREADS \
                --minMappingQuality 40 \
                --region chr10 \
                --binSize 1000 \
                --bamfiles SRR1274742_bwameth_default.bam SRR1274743_bwameth_default.bam \
                SRR1274744_bwameth_default.bam SRR1274745_bwameth_default.bam \
                SRR1653162_bwameth_default.bam SRR2878513_bwameth_default.bam\
                SRR2878520_bwameth_default.bam \
                --labels  SRR1274742_bw SRR1274743_bw \
                SRR1274744_bw SRR1274745_bw \
                SRR1653162_bw SRR2878513_bw\
                SRR2878520_bw \
                -out test_readCounts.npz \
                --outRawCounts test_readCounts.tab

plotCorrelation \
    -in test_readCounts.npz \
    --corMethod spearman --skipZeros \
    --removeOutliers \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o test_heatmap_spearmancorr_readCounts.png   \
    --outFileCorMatrix test_spearmancorr_readcounts.tab

plotPCA -in test_readCounts.npz \
        -o test_PCA_readCounts.png \
        -T "PCA of read counts"

# let's include bismarkk stuff as well
# to be able to fetch data from baubec let-s fetch data from the baubeclabserver

# mkdir -p ~/cg_context_baubec

# scp -i ~/.ssh/cloudServer.key ubuntu@172.23.28.228:/home/ubuntu/MOUNT2/imallona/cg_context/*bam .


# but rather do it with the real datasets
# do they map differently?

# # paired end
# 20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40.bam
# 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam


# ## beware not map40 here?
# SRR1653150_1_cutadapt_sickle_bismark_bt2.deduplicated.bam
# SRR1653150_1_cutadapt_sickle_bwameth_default_dup_marked.bam


for fn in  "20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40.bam" \
               20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam \
               SRR1653150_1_cutadapt_sickle_bismark_bt2.deduplicated.bam \
               SRR1653150_1_cutadapt_sickle_bwameth_default_dup_marked.bam
do
    
    echo $fn
    if [ ! -f "$fn".bai ]
    then
        samtools sort --threads $NTHREADS -o "$fn".sorted $fn
        mv "$fn".sorted $fn
        samtools index $fn
    fi
done

               
multiBamSummary bins \
                --numberOfProcessors $NTHREADS \
                --minMappingQuality 40 \
                --region chr10 \
                --binSize 1000 \
                --bamfiles \
                "20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40.bam" \
                20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam \
                SRR1653150_1_cutadapt_sickle_bismark_bt2.deduplicated.bam \
                SRR1653150_1_cutadapt_sickle_bwameth_default_dup_marked.bam \
                --labels \
                20151223.B-MmES_TKOD3A1c1_bt2 \
                20151223.B-MmES_TKOD3A1c1_bw \
                SRR1653150_bt2 \
                SRR1653150_bw \
                -out test_readCounts.npz \
                --outRawCounts test_readCounts.tab

plotCorrelation \
    -in test_readCounts.npz \
    --corMethod spearman --skipZeros \
    --removeOutliers \
    --plotTitle "Spearman Correlation of Read Counts" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o test_heatmap_spearmancorr_readCounts.png   \
    --outFileCorMatrix test_spearmancorr_readcounts.tab

plotPCA -in test_readCounts.npz \
        -o test_PCA_readCounts.png \
        -T "PCA of read counts"



# do they methylate differently?

# awk the meth calls and correlate them
# /home/imallona/mnt/baubec/cg_context/20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked_stranded.txt
# /home/imallona/mnt/baubec/cg_context/SRR2878513__bwameth_default_stranded.txt
# /home/imallona/mnt/baubec/cg_context/SRR2878513__bwameth_default_stranded.txt


deactivate
