#!/bin/bash
##
## methyldackel calls for non cpg methylation
## outputs are processed by cph_calls.Rmd
##
# Izaskun Mallona
# aug the 22nd
# GPL

export HOME=/home/imallona
export TASK="cg_context"
export WD="$HOME"/mnt/nfs/"$TASK"
export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=32

export MAPQ_THRES=40
export MIN_DEPTH=10

export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

cd $WD

## samples to call the non cpg methylation from
## here rather than a report we require a minimum depth etc from the very beginning

cat << EOF  > for_cph.conf
20151223.B-MmES_TKOD3A1c1-3_R,tko+d3a1,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_001,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_001,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_002,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_003,tko+d3a2,bwa_hiseq2k
SRR1274743,tko+d3a2,bwa_miseq_pe
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_001,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_001,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_002,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_003,tko+d3b1,bwa_hiseq2k
SRR1274744,tko+d3b1,bwa_miseq_pe
SRR1274745,tko+d3b1,bwa_miseq_pe
SRR1653150,qko+d3b1,bwa_hiseq2k
SRR1653151,qko+d3b1,bwa_hiseq2k
SRR1653152,qko+d3b1,bwa_hiseq2k
SRR1653153,qko+d3b1,bwa_hiseq2k
SRR1653154,qko+d3b1,bwa_hiseq2k
SRR1653155,qko+d3b1,bwa_hiseq2k
SRR1653156,qko+d3b1,bwa_hiseq2k
SRR1653157,qko+d3b1,bwa_hiseq2k
SRR1653158,qko+d3b1,bwa_hiseq2k
SRR1653159,qko+d3b1,bwa_hiseq2k
SRR1653160,qko+d3b1,bwa_hiseq2k
SRR1653161,qko+d3b1,bwa_hiseq2k
SRR1653162,qko+d3b1,bwa_miseq_pe 
EOF



# while IFS='' read -r line || [[ -n "$line" ]]
# do
#     bam=""
#     cd "$WD"

#     sample=$(echo $line | cut -f1 -d ',')

#     bam=$(find $WD -name "*""$sample""*bwameth_default.bam")
    
#     echo $sample
    
#     if [ -f $bam ]
#     then
#         echo 'methyldackel call here'
#     else
#         echo 'bam not found for ' "$sample"
#     fi
    
# done < for_cph.conf


for genotype in $(cut -f2 -d"," for_cph.conf | sort | uniq)
do
    echo "$genotype"
    current=$(fgrep "$genotype" for_cph.conf | cut -f1 -d"," | paste -d" " -s)
    # currarray=($current)
    # for sample in "${currarray[@]}"
    rm -f tmp
    for sample in $current
    do
        echo $sample
        # finding the bamfiles to be merged afterwards
        bam=$(find $WD -name "*""$sample""*bwameth_default.bam")
        echo $bam >> tmp        
    done

    mkdir -p merged_"$genotype"
    
    samtools merge -@ $NTHREADS -b tmp merged_"$genotype"/"$genotype"_merged.bam

    ## mapq over 40
    samtools view -@ $NTHREADS -h -b -q "$MAPQ_THRES" merged_"$genotype"/"$genotype"_merged.bam \
             -o mapq.bam
    mv -f mapq.bam merged_"$genotype"/"$genotype"_merged.bam

    mv -f tmp merged_"$genotype"/"$genotype"_merged.bam.sources

    # methyldackel call

    samtools index -@ $NTHREADS \
             merged_"$genotype"/"$genotype"_merged.bam \
             merged_"$genotype"/"$genotype"_merged.bam.bai

    report=merged_"$genotype"/"$genotype"_ch_d_"$MIN_DEPTH"_mapq_"$MAPQ_THRES".cytosine_report.txt

    # this call does not take into account the min depth nor the mapq thres! methyldackel
    # does not (might not?) filter for that without the mergecontext stuff
    $METHYLDACKEL extract \
                  -q $MAPQ_THRES \
                  -@ $NTHREADS \
                  -d $MIN_DEPTH \
                  --cytosine_report \
                  --CHH \
                  --CHG \
                  -o merged_"$genotype"/"$genotype"_ch_d_"$MIN_DEPTH"_mapq_"$MAPQ_THRES" \
                  $MM9 \
                  merged_"$genotype"/"$genotype"_merged.bam

    # this only for chromosome 1 for testing purposes (to be parallelized for all
    # chromosomes @todo)
    # grep -w "^chr1" merged_"$genotype"/"$genotype"_ch_d_"$MIN_DEPTH"_mapq_"$MAPQ_THRES".cytosine_report.txt > tmpchr1
    # mv -f tmpchr1 merged_"$genotype"/"$genotype"_ch_d_"$MIN_DEPTH"_mapq_"$MAPQ_THRES".cytosine_report.txt

    # sampled 10000 rows
    # head -1000000 "$report" > tmphead
    # mv -f tmphead "$report"
    
    
#     # only stuff that is methylated at least once and with a coverage of at least 5
#     echo 'getting 10k covered Cs only'
#     awk '{
# FS=OFS="\t"; 
# if ($4 >= 1 && $4+$5 >= 5)
#  print $0
# }' \
#         "$report" | head -100000  > tmp_covered

    
    # only stuff  with a coverage of at least 10
    # echo 'getting 10k covered Cs only'
    awk '{
FS=OFS="\t"; 
if ($4+$5 >= 10)
 print $0
}' \
        "$report" > tmp_covered

     mv -f tmp_covered "$report"
     
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
          "$report"_cytosine_report_slop.fa > tmp

    rm -rf "$report"_cytosine_report_slop.fa
    # get meth and unmeth stats
    
    # split into cg and non cg 
    awk '{OFS=FS="\t"; if ($6 == "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' tmp > cg
    awk '{OFS=FS="\t"; if ($6 != "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' tmp > ch
    
    # split into meth and unmeth
    awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' cg > cg_meth
    awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' cg > cg_unmeth
        
    awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' ch > ch_meth
    awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' ch > ch_unmeth
    wc -l cg_meth cg_unmeth ch_meth ch_unmeth

    # count instancees by motif
    for item in  cg_meth cg_unmeth ch_meth ch_unmeth
    do
        cut -f9 $item | sort | uniq -c | sed 's/^ *//' > \
                                             merged_"$genotype"/motif_counts_"$item".txt
    done

    rm -f cg ch cg_meth cg_unmeth ch_meth ch_unmeth

    ## this naming is inconsitent, @todo fix
    mv tmp  merged_"$genotype"/raw_report_"$item".txt
    gzip  merged_"$genotype"/raw_report_"$item".txt



done


# ## adding extra summaries with different depths
# ## i.e. higher coverages than 5 reads
# ## @todo code this; replaced the upper chunk to go for 10 reads


# for genotype in $(cut -f2 -d"," for_cph.conf | sort | uniq)
# do
#     cd $WD
#     echo "$genotype"
#     current=$(fgrep "$genotype" for_cph.conf | cut -f1 -d"," | paste -d" " -s)
#     # currarray=($current)

#     item="ch" # this naming is inconsistent, @todo fix
#     for depth in 5 10 25 50
#     do
#         cd $WD/merged_"$genotype"
#         depth_path=merged_"$genotype"/min_depth_"$depth"
#         mkdir -p "$deph_path"
#         cd "$depth_path"
        
#         zcat ../merged_"$genotype"/raw_report_"$item".txt.gz | \
#             awk -v depth="$depth" '{
# FS=OFS="\t"; 
# if ($4+$5 >= depth)
#  print $0
# }' > tmp

#         # # split into cg and non cg 
#         # awk '{OFS=FS="\t"; if ($6 == "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' tmp > cg
#         # awk '{OFS=FS="\t"; if ($6 != "CG") print $1,$2,$3,$4,$5,$6,$7,$8,toupper($9)}' tmp > ch
        
#         # # split into meth and unmeth
#         # awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' cg > cg_meth
#         # awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' cg > cg_unmeth
        
#         # awk '{OFS=FS="\t"; if ($4 > 0) print $0 }' ch > ch_meth
#         # awk '{OFS=FS="\t"; if ($4 == 0) print $0 }' ch > ch_unmeth
#         # wc -l cg_meth cg_unmeth ch_meth ch_unmeth

#         # # count instancees by motif
#         # for item in  cg_meth cg_unmeth ch_meth ch_unmeth
#         # do
#         #     cut -f9 $item | sort | uniq -c | sed 's/^ *//' > \
#         #                                          "$depth_path"/motif_counts_"$item".txt
#         # done

#         # rm -f cg ch cg_meth cg_unmeth ch_meth ch_unmeth
#         cd $WD
#     done
# done
