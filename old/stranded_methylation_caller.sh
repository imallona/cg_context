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
NTHREADS=20

QUALIMAP="$SOFT"/qualimap/qualimap_v2.2.1/qualimap
METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
MARKDUPLICATES=/usr/local/software/picard-tools-1.96/MarkDuplicates.jar
BEDTOOLS="$SOFT"/bedtools/bin/bedtools

MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa

cd $WD

## even though this is (probably) not well mapped since 50 nt long and run with bwa-meth bwa-mem default
toy=SRR1653150_1_cutadapt_sickle_bwameth_default.bam
samtools sort $toy --threads $NTHREADS > sorted.bam
mv sorted.bam $toy

"$QUALIMAP"  bamqc \
             -bam "$toy" \
             -gd mm9 \
             -outdir "$(basename $toy .bam)"_qualimap \
             --java-mem-size=10G \
             -nt $NTHREADS
            

## mark duplicates

## beware of the ram
java -jar -XX:ParallelGCThreads=$NTHREADS \
     $MARKDUPLICATES INPUT=$WD/"$toy" \
     OUTPUT=$WD/"$(basename $toy .bam)""_dup_marked.bam" \
     METRICS_FILE=$WD/"$(basename $toy .bam)""_dup_marked.metrics"

## same for the others

for toy in 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default.bam SRR2878513__bwameth_default.bam
do
    samtools sort $toy --threads $NTHREADS > sorted.bam
    mv sorted.bam $toy
    
    "$QUALIMAP"  bamqc \
             -bam "$toy" \
             -gd mm9 \
             -outdir "$(basename $toy .bam)"_qualimap \
             --java-mem-size=40G \
             -nt $NTHREADS


    java -jar -XX:ParallelGCThreads=$NTHREADS \
         $MARKDUPLICATES INPUT=$WD/"$toy" \
         OUTPUT=$WD/"$(basename $toy .bam)""_dup_marked.bam" \
         METRICS_FILE=$WD/"$(basename $toy .bam)""_dup_marked.metrics"
    
done


# split for/rev strands

## this might not be easy for paired end stuff
# https://broadinstitute.github.io/picard/explaian-flags.html
# http://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/

## methyldackel here
\
# MAPQ_THRES=40
# this to extract watson ones
# samtools view -f 99 -f 147 -bq $MAPQ_THRES file.bam > watson.bam


# algorithm start
#   retrieve 99 and 147 (watsons)

#   retrieve 83 and 163 (cricks)

#   for strand in watson or crick
#       call methyldackel
#       check for variants
      
#    for end
          
# algorithm end

# SRR1653150_1_cutadapt_sickle_bwameth_default_dup_marked.bam 
# for bam in 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam \
#                SRR2878513__bwameth_default_dup_marked.bam


echo 'paired end stuff'

MAPQ_THRES=40
NTHREADS=10

for bam in 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam
do

    echo $bam
    fn="$WD"/$(basename $bam .bam)
    mkdir -p $fn
    
    # proper pairs
    # samtools view -@ $NTHREADS -f 163 -f 83 -bq $MAPQ_THRES $bam > "$fn"/crick_paired.bam
    # samtools view -@ $NTHREADS -f 99 -f 147 -bq $MAPQ_THRES $bam > "$fn"/watson_paired.bam

    samtools view -@ $NTHREADS -f 163  -bq $MAPQ_THRES $bam > "$fn"/crick_paired_163.bam
    samtools view -@ $NTHREADS -f 83 -bq $MAPQ_THRES $bam > "$fn"/crick_paired_83.bam
    
    samtools view -@ $NTHREADS -f 99 -bq $MAPQ_THRES $bam > "$fn"/watson_paired_99.bam
    samtools view -@ $NTHREADS -f 147 -bq $MAPQ_THRES $bam > "$fn"/watson_paired_147.bam
    
    
    # http://www.samformat.info/sam-format-flag
    # crick 89 153
    samtools view -@ $NTHREADS -f 89 -bq $MAPQ_THRES $bam > "$fn"/crick_89.bam
    samtools view -@ $NTHREADS -f 153 -bq $MAPQ_THRES $bam > "$fn"/crick_153.bam
    
    # watson 165 101
    samtools view -@ $NTHREADS -f 169 -bq $MAPQ_THRES $bam > "$fn"/watson_169.bam
    samtools view -@ $NTHREADS -f 105 -bq $MAPQ_THRES $bam > "$fn"/watson_105.bam
    
done

echo 'on methyldackel'

for bam in 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam
do
    echo $bam
    fn="$WD"/$(basename $bam .bam)

    for subset in $(find $fn -name "*bam")
    do
        echo $subset
        cd $fn
           
        $METHYLDACKEL mbias  --txt -@ $NTHREADS $MM9 $subset $(basename $subset .bam) 
                      
        $METHYLDACKEL extract \
                      -@ $NTHREADS \
                      -d 5 \
                      -D 2000 \
                      --CHG \
                      --CHH \
                      $MM9 \
                      $subset \
                      -o $(basename $subset .bam) 

        cd ..
        
    done
done


# out of curiosity, overlaps

head -10000 watson_paired_147_CpG.bedGraph > watson147.head
head -10000 crick_paired_83_CpG.bedGraph > crick83.head

bedtools intersect -sorted \
         -a watson147.head \
         -b crick83.head | wc -l

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome

bedtools slop -i  watson147.head -g mm9.genome -l 0 -r 0 | 
    bedtools intersect \
             -wa -wb \
             -sorted \
             -a - \
             -b crick83.head | head


bedtools slop -i  watson147.head -g mm9.genome -l 0 -r 1 | 
    bedtools intersect \
             -wa -wb \
             -sorted \
             -a - \
             -b crick83.head | head -2

fgrep 3008937 watson147.head

## ok, makes sense

## let's extract the same for each C


for bam in 20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked.bam
do
    echo $bam
    fn="$WD"/$(basename $bam .bam)

    for subset in $(find $fn -name "*bam")
    do
        echo $subset
        cd $fn

        $METHYLDACKEL extract \
                      -@ $NTHREADS \
                      --cytosine_report \
                      $MM9 \
                      $subset \
                      -o $(basename $subset .bam) 

        cd ..
        
    done
done

## actually this is all stupid, better do that directly to the first bam file,
## without flag parsing

# coding the simple_stranded_dna_meth_caller instead
