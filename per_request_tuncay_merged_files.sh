#!/bin/bash
##
## Per request Baubec
# Hi Izaskun, do you have a file where you have the M and U counts for each CpG position
# in the genome,
# including surrounding NNNCGNNN motif and strand info?
# For the merged TKO+d3a2 and TKO+d3b1 samples.
# Chr19-only would be sufficient if that file would be too large to work with.
# I want to test a few things regarding how to calculate preferences according to coverages, etc…
# and if it would be feasible to work on individual sites, rather than average across the genome.
##
## 03 dec 2019
## Izaskun Mallona
## GPL v3 license

export HOME=/home/imallona
export TASK="cg_context"
export SUBTASK="per_tuncay_request_no_coverage_filtering"

export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=20
export MAPQ_THRES=40
export MIN_DEPTH=1 # not used, but anyway

export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export WD="$HOME"/"$TASK"/"$SUBTASK"
export CONFIG_FILE="$TASK".conf

process(){
    bam="$1"
    nthreads="$2"
    mincov="$3"
    BEDTOOLS="$4"
    METHYLDACKEL="$5"

    bam=$(basename $bam)
    sample=$(basename $bam .bam)
    samtools index -@ $nthreads "$(basename $bam)" "$(basename $bam)".bai

    ## mind the mapq40 filtering was done before
    $METHYLDACKEL extract \
                  -@ $nthreads \
                  --cytosine_report \
                  --CHH \
                  --CHG \
                  -o "$(basename $bam)"_ch_d \
                  $MM9 \
                  "$(basename $bam)"

    report="$bam"_ch_d.cytosine_report.txt

    pigz --processes $nthreads $report
    
    ## mincov
    zcat "$report".gz | awk -v mincov="$mincov" '{
FS=OFS="\t"; 
if ($4+$5 >= mincov)
 print $0
}' | gzip > tmp_covered_"$sample"

    mv -f tmp_covered_"$sample" "$report".gz

    ## motif retrieval
    zcat "$report".gz | awk '
{ 
  OFS=FS="\t";
   print $1,$2-1,$2,$4,$5,$3,$7;
}
'  |
        "$BEDTOOLS" slop -i - \
                    -g mm9.genome \
                    -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$report"_cytosine_report_slop.fa \
                    -tab \
                    -s
    
    zcat "$report".gz | paste - "$report"_cytosine_report_slop.fa > tmp_"$sample"

    fgrep chr19 tmp_"$sample" > chr19_"$sample"_tuncay.txt
    mv tmp_"$sample" whole_"$sample"_tuncay.txt

    pigz --processes $nthreads chr19_"$sample"_tuncay.txt
    pigz --processes $nthreads whole_"$sample"_tuncay.txt
 
}


mkdir -p $WD
cd $WD

## these come from recap_by_motif_frequencies.sh
cat << EOF  > "$CONFIG_FILE"
tko+d3a2_merged,tko+d3a2,merged
tko+d3b1_merged,tko+d3b1,merged
EOF

# bam files location
ln -s ~/mnt/nfs/cg_context

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
  
    bam=$(find -L $WD -type f -name "*bam" | fgrep "$sample"| head -1)

    echo $bam
    
    mkdir "$sample"
    cd "$_"
    
    ln -s $bam

    mysql --user=genome \
      --host=genome-euro-mysql.soe.ucsc.edu -A -P 3306 \
      -e "select chrom, size from mm9.chromInfo" > mm9.genome

    process $(basename $bam) \
         $NTHREADS \
         1 \
         $BEDTOOLS \
         $METHYLDACKEL | tee -a no_filtering.log

    cd ..
    
done < "$CONFIG_FILE"

