#!/bin/bash
##
## methyldackel calls for non cpg methylation
## outputs are processed by cph_quickchanges_simpler.Rmd (i.e. bamfiles split by sequencing)
##
# Izaskun Mallona
## 23rd oct
# GPL

export HOME=/home/imallona
export TASK="cg_context"
export WD="$HOME"/mnt/nfs/"$TASK"/split_control
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=24
export MIN_DEPTH=10
export MAPQ_THRES=40

export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

mkdir -p $WD

cd $WD

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome




## bulk bamfiles check

for fn in $(find .. -name "*bam" ! -name "*dup_included*" ! -name "*with_duplicates*" \
                  ! -name "*merged*")
do
    bam=$(basename $fn)
    sample=$(basename $bam .bam)
    echo $sample    
    ###
    
    echo "$(date) cph calls for sample $sample"
    
    # ## mapq over 40
    # samtools view -@ $NTHREADS -h -b -q "$MAPQ_THRES" "$bam" \
    #          -o "$sample"_mapq.bam
    
    # mv -f "$sample"_mapq.bam "$bam"

    ln -s $fn
    
    samtools index -@ $NTHREADS "$bam" "$bam".bai

    # this call does not take into account the min depth nor the mapq thres! methyldackel
    # does not (might not?) filter for that without the mergecontext flag
    $METHYLDACKEL extract \
                  -q $MAPQ_THRES \
                  -@ $NTHREADS \
                  -d $MIN_DEPTH \
                  --cytosine_report \
                  --CHH \
                  --CHG \
                  -o "$bam"_ch_d_"$MIN_DEPTH"_mapq_"$MAPQ_THRES" \
                  $MM9 \
                  "$bam"

    report="$bam"_ch_d_"$MIN_DEPTH"_mapq_"$MAPQ_THRES".cytosine_report.txt

    # only stuff  with a coverage of at least 10
    # echo 'getting 10k covered Cs only'
    awk '{
FS=OFS="\t"; 
if ($4+$5 >= 10)
 print $0
}' \
        "$report" > tmp_covered_"$sample"

    mv -f tmp_covered_"$sample" "$report"
    
    ## motif retrieval

    awk '
{ 
  OFS=FS="\t";
   print $1,$2-1,$2,$4,$5,$3,$7;
}
' "$report"  |
        "$BEDTOOLS" slop -i - \
                    -g "$WD"/mm9.genome \
                    -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$report"_cytosine_report_slop.fa \
                    -tab \
                    -s
    paste "$report" \
          "$report"_cytosine_report_slop.fa > tmp_"$sample"

    rm -rf "$report"_cytosine_report_slop.fa
    # get meth and unmeth stats
    
    # split into cg and non cg 
    awk '{OFS=FS="\t"; if ($6 == "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' \
        tmp_"$sample" > cg_"$sample"
    awk '{OFS=FS="\t"; if ($6 != "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' \
        tmp_"$sample" > ch_"$sample"
    
    # split into meth and unmeth
    awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' cg_"$sample" > cg_meth_"$sample"
    awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' cg_"$sample" > cg_unmeth_"$sample"
    
    awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' ch_"$sample" > ch_meth_"$sample"
    awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' ch_"$sample" > ch_unmeth_"$sample"
    wc -l cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"

    # count instancees by motif
    for item in  cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"
    do
        cut -f9 "$item" | sort | uniq -c | sed 's/^ *//' > \
                                               "$sample"_motif_counts_"$item".txt
    done

    rm -f cg_"$sample" ch_"$sample" cg_meth_"$sample" cg_unmeth_"$sample" \
       ch_meth_"$sample" ch_unmeth_"$sample"

    mv tmp_"$sample" "$sample"_raw_report.txt
    rm -f "$sample"_raw_report.txt

done 

## cleaning (removal of .bai s etc)
