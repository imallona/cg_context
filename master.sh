#!/bin/bash

# preliminar.sh                       contains qualcheck etc using bwameth
# bowtie_based.sh                     same with bowtie2
# stranded_dna_meth_caller.sh         meth calls using methyldackel using really convoluted sam flags parsing
# simple_stranded_dna_meth_caller.sh  kiss version of the upper
# simple_stranded_dna_meth_caller.R   analyzes the upper
# paired_end_workflow                 detected that a sample did not map as pe, writing the whole flow


# this is getting messy, let's clarify
# bowtie_based.sh                         bismark stuff
# paired_end_workflow_bulk.sh             bwa meth paired end stuff
# deptools_comparison.sh                  test on deptools representations just in case there were
# bwa and bowtie mapping differences
# diagnostic_plots.R                      preliminar integration
# new_tuncay_paired_data.sh               on dnmts data

# noticed on may the 4th that the slopbed stuff crashes (also does in beech) sometimes, randomly,
# beware; need to write a homemade parser
# apparently not, it works now (??)

# integration.R                           as if diagnostic_plots.R

# stadler
# stadler_single_end.sh                  stadlers'

# set up an NFS, let's sync the data there
for fn in $(find ~/cg_context -name "*stranded.txt.gz")
do
    rsync -avt $fn ~/mnt/nfs/cg_context/ &
done

for fn in $(find ~/cg_context_bulk -name "*stranded.txt.gz")
do
    rsync -avt $fn ~/mnt/nfs/cg_context/ &
done

for fn in $(find ~/cg_context_new_tuncay/ -name "*stranded.txt.gz")
do
    rsync -avt $fn ~/mnt/nfs/cg_context/ &
done

## even DMMD NAS and S3IT NFS sync (not needed)
for fn in $(find ~/mnt/baubec/imallona2mmanzo/ -name "*stranded.txt.gz")
do
    #echo $fn
    rsync -avt $fn ~/mnt/nfs/cg_context/ &
done

chmod 644 ~/mnt/nfs/cg_context/*gz

# consistency checks to count lines

for fn in $(find ~/mnt/nfs/cg_context -name "*strand**txt.gz")
do
    nlines=$(zcat $fn | wc -l | cut -f1 -d" ")
    
    if [[ ! $nlines -eq 21722957 ]] 
    then 
        echo 'consistency check failed for ' "$fn"
        echo $nlines
    fi
done

# crashes for two
# ^[[Aconsistency check failed for  /home/imallona/mnt/nfs/cg_context/SRR299053_bwameth_default_stranded.txt.gz
# 23729652
# ^[[Bconsistency check failed for  /home/imallona/mnt/nfs/cg_context/SRR299054_bwameth_default_stranded.txt.gz
# 23733897


export HOME=/home/imallona
export TASK="cg_context"
export WD="$HOME"/"$TASK"
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa

export MAPQ_THRES=40
export NTHREADS=6
export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools


cd $WD
for sample in SRR299053 SRR299054
do
    cd $sample
    
    bam="$sample"_bwameth_default.bam

    $METHYLDACKEL extract \
                  -q $MAPQ_THRES \
                  -@ $NTHREADS \
                  --cytosine_report \
                  $MM9 \
                  $bam \
                  -o $(basename $bam .bam)

    echo 'the bedtools slop migth be broken'
    
    awk '
{ 
  OFS=FS="\t";
   print $1,$2-1,$2,$4,$5,$3,$7;
}
' "$(basename $bam .bam)".cytosine_report.txt  |
        "$BEDTOOLS" slop -i - \
                    -g "$WD"/mm9.genome \
                    -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$(basename $bam .bam)"_cytosine_report_slop.fa \
                    -tab \
                    -s
    paste "$(basename $bam .bam)".cytosine_report.txt \
          "$(basename $bam .bam)"_cytosine_report_slop.fa > tmp

    rm -rf "$(basename $bam .bam)"_cytosine_report_slop.fa
    ## now get odd and even lines
    awk '{printf "%s%s",$0,(NR%2?FS:RS)}' tmp > bar
    rm -f tmp 
    mv -f bar "$(basename $bam .bam)"_stranded.txt

    echo "$(date) Processing sample $sample ended"

    gzip "$(basename $bam .bam)"_stranded.txt

    rsync -avt "$(basename $bam .bam)"_stranded.txt.gz ~/mnt/nfs/cg_context/

    cd "$WD"
        
done
