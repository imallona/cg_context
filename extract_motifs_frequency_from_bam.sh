#!/bin/bash
##
## Generates discretized DNA methylation count tables from bwameth-aligned bam files
##
## Instances are loci (CpG or CpH) and their surrounding sequence up to a 8-mer motif.
##   E.g. AAACAGGG centered in chrom1:235
##
## Aggregated count tables depict the number of instances that share (in a strand-specific way)
##   the motif, e.g. two instances of AAACAGGG, one centered in chrom1:235 and another in chrX:10
##
## DNA methylation is reported discretizing the methylation value.
##  That is, reports how many DNA methylation instances are fully unmethylated
##  in a bam file, how many have less than 10% methylated reads, how many from 10 to 20% and so on.
##
## Deprecated; extract_motifs_frequency_from_bam_binary.sh should be used instead.
##
## 17th Dec 2018
##
## Izaskun Mallona
## GPLv2

usage(){
    echo ""
    echo "Usage: bash ""$(basename "$0")"" -b bamfile -t nthreads"
    echo ""
    echo "Params"
    echo " -b --bam        WGBS bamfile          [mandatory]"
    echo " -t --threads    number of cores       [defaults to 4]"
    echo " --bedtools      path to bedtools      [defaults to bedtools]"
    echo " --methyldackel  path to methyldackel  [defaults to methyldackel]"
    echo ""
}

# @todo split into items! parallelize by chromosome!
process(){
    bam="$1"
    nthreads="$2"
    BEDTOOLS="$3"
    METHYLDACKEL="$4"
    
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

    ## min coverage 10, @todo better not hardcoded
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

#     awk '{
# OFS=FS="\t"; 
# if ( ($6 == "CG") && (($4/($4+$5)) < 0.2) )  
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_02"
# else if ( ($6 == "CG") && (($4/($4+$5)) < 0.4) ) 
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_04"
# else if ( ($6 == "CG") && (($4/($4+$5)) < 0.6) )
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_06"
# else if ( ($6 == "CG") && (($4/($4+$5)) < 0.8) ) 
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_08"
# else if ( ($6 == "CG") && (($4/($4+$5)) <= 1) ) 
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_1"
# else if ( ($6 != "CG") && (($4/($4+$5)) < 0.2) ) 
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_02"
# else if ( ($6 != "CG") && (($4/($4+$5)) < 0.4) ) 
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_04"
# else if ( ($6 != "CG") && (($4/($4+$5)) < 0.6) )
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_06"
# else if ( ($6 != "CG") && (($4/($4+$5)) < 0.8) ) 
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_08"
# else if ( ($6 != "CG") && (($4/($4+$5)) <= 1) ) 
#   print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_1"
# else 
#   print $0,$4/($4+$5) > "discretized/error"
# }' tmp_"$sample"

        awk '{
OFS=FS="\t"; 
if ( ($6 == "CG") && ($4 == 0) )  
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_000"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.1) ) 
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_010"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.2) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_020"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.3) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_030"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.4) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_040"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.5) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_050"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.6) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_060"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.7) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_070"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.8) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_080"
else if ( ($6 == "CG") && (($4/($4+$5)) < 0.9) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_090"
else if ( ($6 == "CG") && (($4/($4+$5)) <= 1) ) 
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/cg_100"

else if ( ($6  != "CG") && ($4 == 0) )  
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_000"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.1) ) 
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_010"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.2) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_020"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.3) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_030"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.4) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_040"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.5) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_050"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.6) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_060"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.7) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_070"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.8) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_080"
else if ( ($6  != "CG") && (($4/($4+$5)) < 0.9) )
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_090"
else if ( ($6  != "CG") && (($4/($4+$5)) <= 1) ) 
  print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9) > "discretized/ch_100"

else 
  print $0,$4/($4+$5) > "discretized/error"
}' tmp_"$sample" 

    ## some cg or ch files might be missing

    for fn in $(find discretized -type f)
    do
        echo $fn
        item=$(basename $fn)
        cut -f9 "$fn" | sort | uniq -c | sed 's/^ *//' > \
                                             discretized/"$sample"_motif_counts_"$item".txt
    done

    rm -f discretized/c*
    mv tmp_"$sample" "$sample"_raw_report.txt
    rm -f "$sample"_raw_report.txt

}

if [[ $# -eq 0 ]] ; then
    usage
    exit 0
fi

bam=
nthreads=4
BEDTOOLS=bedtools
METHYLDACKEL=methyldackel
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
         --bedtools)           shift
                               BEDTOOLS=$1
                               ;;
         --methyldackel)       shift
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


echo "$(date)" Processing sample $sample from "$bam" with "$nthreads" threads start 

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" > mm9.genome

# # bam=test.bam
# bash extract_motifs_frequency_from_bam.sh -b test.bam \
#      -t 4 \
#      -s $(basename test.bam .bam)
#      --bedtools $BEDTOOLS
process $bam $nthreads $BEDTOOLS $METHYLDACKEL

echo "$(date)" Processing sample $sample from "$bam" with "$nthreads" threads end
