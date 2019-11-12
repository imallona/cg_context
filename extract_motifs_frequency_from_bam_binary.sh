#!/bin/bash
##
## 04 nov 2019
## Extracts meth and un meth motifs from bisulfite bam files, and counts them
## 
## Assumes is mm9! hardcoded
##
## Izaskun Mallona

MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa

usage(){
    echo ""
    echo "Usage: bash ""$(basename "$0")"" -b bamfile -t nthreads"
    echo ""
    echo "Params"
    echo " -b --bam        WGBS bamfile          [mandatory]"
    echo " -t --threads    number of cores       [defaults to 4]"
    echo " -c --coverage   min num reads         [defaults to 10 nonstrandspecific]"
    echo " --bedtools      path to bedtools      [defaults to bedtools]"
    echo " --methyldackel  path to methyldackel  [defaults to MethylDackel]"
    echo ""
}

# @todo split into items! parallelize by chromosome!
process(){
    bam="$1"
    nthreads="$2"
    mincov="$3"
    BEDTOOLS="$4"
    METHYLDACKEL="$5"
    
    sample=$(basename $bam .bam)
    samtools index -@ $nthreads "$bam" "$bam".bai

    ## mind the mapq40 filtering was done before
    $METHYLDACKEL extract \
                  -@ $nthreads \
                  --cytosine_report \
                  --CHH \
                  --CHG \
                  -o "$bam"_ch_d \
                  $MM9 \
                  "$bam"

    report="$bam"_ch_d.cytosine_report.txt

    ## mincov
    awk -v mincov="$mincov" '{
FS=OFS="\t"; 
if ($4+$5 >= mincov)
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
                    -g mm9.genome \
                    -l 3 -r 4 -s | \
        "$BEDTOOLS" getfasta -fi $MM9 \
                    -bed - \
                    -fo "$report"_cytosine_report_slop.fa \
                    -tab \
                    -s
    paste "$report" \
          "$report"_cytosine_report_slop.fa > tmp_"$sample"

    rm -rf "$report"_cytosine_report_slop.fa

    ## iterate over diff meth states
    mkdir -p discretized

    
    # split into cg and non cg

    # ## old start
    # awk '{OFS=FS="\t"; if ($6 == "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' \
    #     tmp_"$sample" > cg_"$sample"
    # awk '{OFS=FS="\t"; if ($6 != "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' \
    #     tmp_"$sample" > ch_"$sample"
    
    # # split into meth and unmeth
    # awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' cg_"$sample" > cg_meth_"$sample"
    # awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' cg_"$sample" > cg_unmeth_"$sample"
    
    # awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' ch_"$sample" > ch_meth_"$sample"
    # awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' ch_"$sample" > ch_unmeth_"$sample"

    # ## old end  
    
    ## split into cg and non cg, meth and unmeth in one run
    awk -v cgmeth=cg_meth_"$sample" -v chmeth=ch_meth_"$sample" \
        -v cgunmeth=cg_unmeth_"$sample" -v chunmeth=ch_unmeth_"$sample" '{
OFS=FS="\t"; 
if ($6 == "CG" && $4 > 0) {
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > cgmeth}
else if ($6 == "CG" && $4 == 0) {
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > cgunmeth}
else if ($6 != "CG" && $4 > 0) {
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > chmeth}
else if ($6 != "CG" && $4 == 0) {
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > chunmeth}
else
  print
}' tmp_"$sample"

    wc -l cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"    

    # count instances by motif
    for item in  cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"
    do
        cut -f9 "$item" | sort | uniq -c | sed 's/^ *//' > \
                                               "$sample"_motif_counts_"$item".txt
    done

    rm -f cg_"$sample" ch_"$sample" cg_meth_"$sample" cg_unmeth_"$sample" \
       ch_meth_"$sample" ch_unmeth_"$sample"

    mv tmp_"$sample" "$sample"_raw_report.txt
    rm -f "$sample"_raw_report.txt

}

if [[ $# -eq 0 ]] ; then
    usage
    exit 0
fi

bam=
nthreads=4
mincov=10
BEDTOOLS=bedtools
METHYLDACKEL=MethylDackel
while [ "$1" != "" ]; do
    case $1 in
        -b | --bam )           shift
                               bam="$1"
                               ;;
        -t | --threads )       shift
                               nthreads=$1
                               ;;
        -s | --sample)         shift
                               sample=$1
                               ;;
        -c | --coverage)       shift
                               mincov=$1
                               ;;
        --bedtools)            shift
                               BEDTOOLS=$1
                               ;;
        --methyldackel)        shift
                               METHYLDACKEL=$1
                               ;;
        -h | --help )          usage
                               exit
                               ;;
        * )                    usage
                               exit 1
    esac
    shift
done

if [[ "$bam" == "" ]] ; then
    usage
    exit 0
fi


echo "$(date)" Processing $sample from "$bam" with "$nthreads" threads and mincov "$mincov" start

mysql --user=genome \
      --host=genome-euro-mysql.soe.ucsc.edu -A -P 3306 \
      -e "select chrom, size from mm9.chromInfo" > mm9.genome

process $bam $nthreads $mincov $BEDTOOLS $METHYLDACKEL

echo "$(date)" Processing $sample from "$bam" with "$nthreads" threads and mincov "$mincov" end
