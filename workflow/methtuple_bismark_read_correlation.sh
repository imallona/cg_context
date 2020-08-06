#!/bin/bash
##
## Evaluates the intra-read CpG - CpH correlation/dependence for one sample (paired end)
##
## Izaskun Mallona
## 5 Aug 2020

##
# for instance SRR1274743 TKO_DNMT3A2_bisulfite_miseq paired end

sample=SRR1274743


export START=/home/Shared_s3it/imallona/cg_context/nar_review
export TASK=autocorrelation
export WD="$START"/"$TASK"

export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=18

export MAPQ_THRES=40

export FASTQC=/usr/local/software/FastQC/fastqc
export SICKLE="$HOME"/soft/sickle/sickle-1.33/sickle
export ##CUTADAPT=/usr/local/bin/cutadapt
export CUTADAPT="$VIRTENVS"/cutadapt/bin/cutadapt
export QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
export MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools
export FASTQDUMP=/usr/local/software/sratoolkit.2.9.0-ubuntu64/bin/fastq-dump # @todo update

export ILLUMINA_UNIVERSAL="AGATCGGAAGAG"
export ILLUMINA="CGGTTCAGCAGGAATGCCGAGATCGGAAGAGCGGTT"

export BOWTIE_PATH=~/soft/bowtie/bowtie-1.2.2/
export BISMARK=~/soft/bismark/Bismark_v0.19.1/bismark

mkdir -p $WD

cd $_


echo 'uncomment to download'
$FASTQDUMP -I --gzip --split-files $sample
export sample="$sample"
echo $sample

# for r in 1 2
# do
#     curr=${WD}/${sample}
#     mkdir -p $curr
    
#     $FASTQC ${WD}/"$sample"_"$r".fastq.gz --outdir ${curr} \
#             -t $NTHREADS &> ${curr}/${sample}_fastqc.log

# done

source $VIRTENVS/cutadapt/bin/activate

cutadapt \
    -j $NTHREADS \
    -b $ILLUMINA_UNIVERSAL -b $ILLUMINA \
    -B $ILLUMINA_UNIVERSAL -B $ILLUMINA \
    -o "$WD"/"${sample}"_1_cutadapt.fastq.gz \
    -p "$WD"/"${sample}"_2_cutadapt.fastq.gz \
    "$sample"_1.fastq.gz "$sample"_2.fastq.gz &> "$WD"/"$sample"_cutadapt.log

deactivate


"$SICKLE" pe \
          -f "$WD"/"$sample"_1_cutadapt.fastq.gz \
          -r "$WD"/"$sample"_2_cutadapt.fastq.gz \
          -o "$WD"/"$sample"_1_cutadapt_sickle.fastq.gz \
          -p "$WD"/"$sample"_2_cutadapt_sickle.fastq.gz \
          -t sanger \
          -s "$WD"/"$sample"_cutadapt_sickle_singles.fastq.gz \
          -g &> "$WD"/"$sample"_cutadapt_sickle.log


rm -f "$WD"/"${sample}"_1_cutadapt.fastq.gz "$WD"/"${sample}"_2_cutadapt.fastq.gz

# for r in 1 2
# do
#     curr="$sample"_"$r"_cutadapt_sickle
#     mkdir -p "$WD"/"$curr"
#     $FASTQC "$WD"/${sample}_"$r"_cutadapt_sickle.fastq.gz \
#             --outdir "$WD"/"$curr" \
#             -t $NTHREADS &> ${curr}/${sample}_"$r"_fastqc.log
# done    

# # bismark here


fw="$sample"_1_cutadapt_sickle.fastq.gz
rv="$sample"_2_cutadapt_sickle.fastq.gz

# # ## bwameth doesn't like reads with different read names, treats them as single end
# # zcat $fw | sed 's/\.1 /\ 1 /g' | gzip -c  > "$fw"_ed.gz
# # zcat $rv | sed 's/\.2 /\ 2 /g'  | gzip -c > "$rv"_ed.gz
# # mv "$fw"_ed.gz  "$fw"
# # mv "$rv"_ed.gz "$rv"

## assumes directional library
( "$BISMARK" --bowtie1 \
             --path_to_bowtie $BOWTIE_PATH \
             --gzip \
             --parallel $NTHREADS \
             --genome $(dirname "$MM9") \
             --1 "$fw" \
             --2 "$rv") 2>&1 | tee "$WD"/"$sample"_bismark.log

# Final Alignment report
# ======================
# Sequence pairs analysed in total:       7893867
# Number of paired-end alignments with a unique best hit: 2933414
# Mapping efficiency:     37.2%

# Sequence pairs with no alignments under any condition:  4758398
# Sequence pairs did not map uniquely:    202055
# Sequence pairs which were discarded because genomic sequence could not be extracted:    4

# Number of sequence pairs with unique best (first) alignment came from the bowtie output:
# CT/GA/CT:       1467637 ((converted) top strand)
# GA/CT/CT:       0       (complementary to (converted) top strand)
# GA/CT/GA:       0       (complementary to (converted) bottom strand)
# CT/GA/GA:       1465777 ((converted) bottom strand)


# # let's try nondirectional
# ( "$BISMARK" --bowtie1 \
#              --path_to_bowtie $BOWTIE_PATH \
#              --gzip \
#              --non_directional \
#              --parallel $NTHREADS \
#              --genome $(dirname "$MM9") \
#              --1 "$fw" \
#              --2 "$rv") 2>&1 | tee "$WD"/"$sample"_bismark_non_directional.log

# ## does not get better, stick to the directional run

# # Final Alignment report
# # ======================
# # Sequence pairs analysed in total:       7893867
# # Number of paired-end alignments with a unique best hit: 2936145
# # Mapping efficiency:     37.2%

# # Sequence pairs with no alignments under any condition:  4748733
# # Sequence pairs did not map uniquely:    208989
# # Sequence pairs which were discarded because genomic sequence could not be extracted:    3

# # Number of sequence pairs with unique best (first) alignment came from the bowtie output:
# # CT/GA/CT:       1466882 ((converted) top strand)
# # GA/CT/CT:       1979    (complementary to (converted) top strand)
# # GA/CT/GA:       2022    (complementary to (converted) bottom strand)
# # CT/GA/GA:       1465262 ((converted) bottom strand)


# remove duplicates

echo dedup

( "$(dirname $BISMARK)"/deduplicate_bismark -p \
                     --output_dir $WD \
                     --bam \
                     "$sample"_1_cutadapt_sickle_bismark_pe.bam ) 2>&1 | \
     tee "$WD"/"$sample"_dedup_bismark.log


# ## run methtuple to 4 tuples at most, including cpg and cph stuff

source ~/virtenvs/methtuple/bin/activate

# Using all paired-end reads, even those that have overlapping mates.
#   However, for the overlapping sequence, only read_1 shall be used
#   (as done by bismark_methylation_extractor --no_overlap).
methtuple "$sample"_1_cutadapt_sickle_bismark_pe.bam \
          --methylation-type CG \
          --methylation-type CHG \
          --methylation-type CHH \
          -m 8 \
          --ignore-duplicates \
          --overlap-filter Bismark

## control only CG
methtuple "$sample"_1_cutadapt_sickle_bismark_pe.bam \
          --methylation-type CG \
          -m 8 \
          --ignore-duplicates \
          --overlap-filter Bismark

deactivate
