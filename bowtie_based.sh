#!/bin/bash
##
## WGBS mapping using bowtie
##
# Izaskun Mallona
# 16th apr 2018
##
## to be run at baubeclab
# GPL

TASK="cg_context"
WD=/home/ubuntu/MOUNT2/imallona/"$TASK"
# DATA="$HOME"/data
DATA=/home/ubuntu/MOUNT2/imallona
SOFT="$HOME"/soft
MM9=/home/ubuntu/REPOS/annotation/mm9_genome/
VIRTENVS=~/virtenvs

NTHREADS_BOWTIE=4
NTHREADS_BISMARK=2

# baubeclab soft
BISMARK=/home/ubuntu/Soft/Bismark-0.19.0/bismark

cd $DATA

echo 'rsync -a -v imallona@imlstaupo.uzh.ch:/home/imallona/cg_context/*cutadapt_sickle.fastq.gz'

mkdir -p $WD
cd $WD

## let's index the genome (rather not, use baubeclab's)
echo single end

for sample in SRR1653150_1
do
    echo $sample
    fastq="$DATA"/"$sample"_cutadapt_sickle.fastq.gz
    ( "$BISMARK" --bowtie2 \
               --genome $MM9 \
               --N 0 \
               -D 20 \
               -R 3 \
               -L 20 \
               --gzip \
               --parallel $NTHREADS_BISMARK \
               -p $NTHREADS_BOWTIE \
               --se "$fastq" ) 2>&1 | tee "$WD"/"$sample"_bismark.log
done

# echo paired end


for sample in SRR2878513_ 20151223.B-MmES_TKOD3A1c1-3_R
do
    echo $sample
    r1="$sample"1_cutadapt_sickle.fastq.gz
    r2="$sample"2_cutadapt_sickle.fastq.gz

    echo $sample
    fastq="$DATA"/"$sample"_cutadapt_sickle.fastq.gz
    ( "$BISMARK" --bowtie2 \
                 --genome $MM9 \
                 --N 0 \
                 -D 20 \
                 -R 3 \
                 -L 20 \
                 --gzip \
                 --parallel $NTHREADS_BISMARK \
                 -p $NTHREADS_BOWTIE \
                 --1 "$r1" \
                 --2 "$r2"
    ) 2>&1 | tee "$WD"/"$sample"_bismark.log
done



# sample=SRR1653150_1_cutadapt_sickle
# # bwameth.py --reference "$MM9" \
# #            "$sample".fastq.gz \
# #            --threads $NTHREADS | samtools view -bS - > "$sample"_bwameth_default.bam 2>&1 | \
# #     tee "$sample"_bwameth_default.log

# # bwameth.py --reference "$MM9" \
# #            "$sample".fastq.gz \
# #            --threads $NTHREADS |  samtools view -bS - > "$sample"_bwameth_default.bam 2> test.log

# ( bwameth.py --reference "$MM9" \
#              "$sample".fastq.gz \
#              --threads $NTHREADS |  samtools view -bS - > "$sample"_bwameth_default.bam ) \
#     3>&1 1>&2 2>&3 | tee "$sample"_bwameth_default.log

# # samtools view SRR1653150_1_cutadapt_sickle_bwameth_default.bam | head -1000 | cut -f2 | sort | uniq -c






# NTHREADS=6
# ## paired end stuff
# for sample in SRR2878513_ 20151223.B-MmES_TKOD3A1c1-3_R 
# do
    
#     source $VIRTENVS/cutadapt/bin/activate

#     cutadapt \
#         -j $NTHREADS \
#         -b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
#         -B $ILLUMINA_UNIVERSAL -B $ILLUMINA \
#         -o "$WD"/"${sample}"1_cutadapt.fastq.gz \
#         -p "$WD"/"${sample}"2_cutadapt.fastq.gz \
#         "$DATA"/"$sample"1.fastq.gz "$DATA"/"$sample"2.fastq.gz &> "$WD"/"$sample"_cutadapt.log

#     deactivate

#     "$SICKLE" pe \
#               -f "$WD"/"$sample"1_cutadapt.fastq.gz \
#               -r "$WD"/"$sample"2_cutadapt.fastq.gz \
#               -o "$WD"/"$sample"1_cutadapt_sickle.fastq.gz \
#               -p "$WD"/"$sample"2_cutadapt_sickle.fastq.gz \
#               -t sanger \
#               -s "$WD"/"$sample"_cutadapt_sickle_singles.fastq.gz \
#               -g &> "$WD"/"$sample"_cutadapt_sickle.log

#     for read in 1 2
#     do
#         curr="$sample"_cutadapt_sickle_"$read"
#         mkdir -p "$WD"/$curr
#         $FASTQC "$WD"/${sample}"$read"_cutadapt_sickle.fastq.gz \
#                 --outdir "$WD"/"$curr" \
#                 -t $NTHREADS &> ${curr}/${sample}"$read"_fastqc.log

#     done
#     # done

#     source $VIRTENVS/bwa-meth/bin/activate

#     fw="$sample"1_cutadapt_sickle.fastq.gz
#     rv="$sample"2_cutadapt_sickle.fastq.gz

#     ( bwameth.py --reference "$MM9" \
#                  $fw $rv \
#                  --threads $NTHREADS |  samtools view -bS - > "$sample"_bwameth_default.bam ) \
#             3>&1 1>&2 2>&3 | tee "$sample"_bwameth_default.log

#     deactivate
# done



# ## extra round of fastqc to detect the 10-mer enrichment
# ## on our samples
# sample=SRR2878513
# for r in in 1 2
# do
#     curr="$sample"_cutadapt_sickle_"$read"_10kmers
#     mkdir -p "$WD"/"$curr"

#     "$FASTQC" "$WD"/"$sample"_"$r"_cutadapt.fastq.gz \
#               -k 10 \
#               --outdir ${curr} \
#               -t $NTHREADS &> ${curr}/${sample}_fastqc.log
# done

# ## check mapq distribution to get rid of multimappers!
