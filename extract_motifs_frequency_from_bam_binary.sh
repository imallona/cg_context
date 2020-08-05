#!/bin/bash
##
## Generates DNA methylation count tables from bwameth-aligned bam files depicting how many
##   instances are methylated (e.g. at least a methylated) read/fragment and how many are not
##   (all reads/fragment unmethylated)
##
## Instances are loci (CpG or CpH) and their surrounding sequence up to a 8-mer motif.
##   E.g. AAACAGGG centered in chrom1:235
##
## Aggregated count tables depict the number of instances that share (in a strand-specific way)
##   the motif, e.g. two instances of AAACAGGG, one centered in chrom1:235 and another in chrX:10
##
## DNA methylation is a binary outcome: either methylated, either non methylated.
##
## A minimum coverage can be requested.
##
## A mask (BED file with intervals to be scrutinized) can be added with \
##  the --intervals flag (added 4th August 2020)
##
## 04 nov 2019
## Izaskun Mallona

MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa

usage(){
    echo ""
    echo "Usage: bash ""$(basename "$0")"" -b bamfile -t nthreads"
    echo ""
    echo "Params"
    echo " -b --bam        WGBS bamfile                      [mandatory]"
    echo " -t --threads    number of cores                   [defaults to 4]"
    echo " -c --coverage   min num reads                     [defaults to 10 nonstrandspecific]"
    echo " --bedtools      path to bedtools                  [defaults to bedtools]"
    echo " --methyldackel  path to methyldackel              [defaults to MethylDackel]"
    echo " -i --intervals  BED file with target intervals    [defaults to None]"
    echo ""
}

# @todo split into items! parallelize by chromosome!
process(){
    bam="$1"
    nthreads="$2"
    mincov="$3"
    BEDTOOLS="$4"
    METHYLDACKEL="$5"
    intervals="$6"

    # bam=$(basename $bam)
    sample=$(basename $bam .bam)
    samtools index -@ $nthreads "$bam" "$bam".bai

    ## mind the mapq40 filtering was done before (bamfile generation)
    $METHYLDACKEL extract \
                  -@ $nthreads \
                  --cytosine_report \
                  --CHH \
                  --CHG \
                  -o "$(basename $bam)"_ch_d \
                  $MM9 \
                  "$bam"

    report="$(basename $bam)"_ch_d.cytosine_report.txt

    ## if no bedfile with intervals to focus in, proceed
    ##   else, intersect the report with the interval BEDFILE
    if [[ "$intervals" != "None"  ]] ; then
        echo 'Applying a mask'        

        ## only grepping elements with data (not 0,0 for meth and unmeth reads)
        awk '{FS=OFS="\t"; print $1,$2,$2+1,$3,$4,$5,$6,$7}' "$report" |
            grep -vP '\t0\t0\t' | \
            "$BEDTOOLS" intersect -a - -b "$intervals" -wa > masked
        
        report="$(basename $bam)"_ch_d.cytosine_report_mask_"$(basename $intervals)".txt
        mv -f masked "$report"
        
    fi

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

    pigz --processes $nthreads tmp_"$sample"

    rm -rf "$report"_cytosine_report_slop.fa

    ## split into cg and non cg, meth and unmeth in one run
    zcat tmp_"$sample".gz | awk -v cgmeth=cg_meth_"$sample" -v chmeth=ch_meth_"$sample" \
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
}' 

    wc -l cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"    

    # count instances by motif
    for item in  cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"
    do
        cut -f9 "$item" | sort | uniq -c | sed 's/^ *//' > \
                                               motif_counts_"$item".txt
    done

    rm -f cg_"$sample" ch_"$sample" cg_meth_"$sample" cg_unmeth_"$sample" \
       ch_meth_"$sample" ch_unmeth_"$sample"

    mv -f tmp_"$sample".gz "$sample"_raw_report.txt.gz
    rm -f "$sample"_raw_report.txt.gz
    rm -f "$bam".bai mm9.genome
    
    pigz --processes $nthreads *motif_counts*

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
intervals=None
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
        -i | --intervals)      shift
                               intervals=$1
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

process $bam $nthreads $mincov $BEDTOOLS $METHYLDACKEL "$intervals"

echo "$(date)" Processing $sample from "$bam" with "$nthreads" threads and mincov "$mincov" end
