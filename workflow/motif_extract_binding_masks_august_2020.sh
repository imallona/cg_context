#!/bin/bash
##
##
## 04 August 2020
## Izaskun mallona


## second motif extract without filtering, and with the extract_motifs_frequency_from_bam_binary.sh
##   updated to contain the mincoverage and masks (e.g. regions to focus on)

export HOME=/home/imallona
export TASK="cg_context"

export SOFT="$HOME"/soft
export MM9=/home/Shared/data/annotation/Mouse/mm9/mm9.fa
export VIRTENVS=~/virtenvs

export NTHREADS=20
export MAPQ_THRES=40
export MIN_DEPTH=0

export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
export BEDTOOLS="$SOFT"/bedtools/bin/bedtools

export EXTRACT_MOTIFS_FREQUENCY_FROM_BAM_BINARY=~/src/cg_context/extract_motifs_frequency_from_bam_binary.sh

export WD='/home/Shared_s3it/imallona/cg_context/nar_review'

export CONFIG_FILE=with_masks.conf


mkdir -p $WD

cd $WD

# folder bams should be there
ls -l bams

ln -s $SRC/cg_context/data
ls -l data

## these come from recap_by_motif_frequencies.sh
cat << EOF  > "$CONFIG_FILE"
20151223.B-MmES_TKOD3A1c1-3_R,tko+d3a1,bwa_hiseq2k
tko+d3a1_merged,tko+d3a1,merged
tko+d3a2_merged,tko+d3a2,merged
tko+d3b1_merged,tko+d3b1,merged
qko+d3b1_merged,qko+d3b1,merged
EOF

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
    # bam='"'"$sample""*bam"'"'
    # find ~/mnt/nfs/cg_context -name $sample -name "*bam"
    bam=$(find -L $WD -type f -name "*bam" | fgrep $sample| head -1)
    echo $bam
done < "$CONFIG_FILE" | wc -l

# ok, they are four

## test with a bam and mask only start

# mkdir test
# cd $_

# samtools view -bs 0.005 "../bams/tko+d3b1_merged.bam" > sampled.bam

# bam="sampled.bam"
# intervals="../data/NEW_DNMT3A1_sites.bed"
# nthreads=20
# mincov=1
# export METHYLDACKEL="$SOFT"/methyldackel/MethylDackel/MethylDackel
# export BEDTOOLS="$SOFT"/bedtools/bin/bedtools


# sample=$(basename $bam .bam)
# samtools index -@ $nthreads "$bam" "$bam".bai

# ## mind the mapq40 filtering was done before (bamfile generation)
# $METHYLDACKEL extract \
#               -@ $nthreads \
#               --cytosine_report \
#               --CHH \
#               --CHG \
#               -o "$(basename $bam)"_ch_d \
#               $MM9 \
#               "$bam"


# report="$bam"_ch_d.cytosine_report.txt

# ## if no bedfile with intervals to focus in, proceed
# ##   else, intersect the report with the interval BEDFILE
# if [[ "$intervals" != "None"  ]] ; then
#     echo 'Applying a mask'        

#     ## only grepping elements with data (not 0,0 for meth and unmeth reads)
#     awk '{FS=OFS="\t"; print $1,$2,$2+1,$3,$4,$5,$6,$7}' "$report" |
#         grep -vP '\t0\t0\t' | \
#         "$BEDTOOLS" intersect -a - -b "$intervals" -wa > masked
    
#     report="$bam"_ch_d.cytosine_report_mask_"$(basename $intervals)".txt
#     mv -f masked "$report"
  
# fi



# bash $EXTRACT_MOTIFS_FREQUENCY_FROM_BAM_BINARY \
#              -b  $bam \
#              -t $NTHREADS \
#              -c 1 \
#              --bedtools $BEDTOOLS \
#              --methyldackel $METHYLDACKEL \
#              --intervals "$mask"| tee -a with_mask_"$(basename $mask .bed)".log



## test with a bam and mask only end


cd $WD

while IFS='' read -r line || [[ -n "$line" ]]
do
    sample=$(echo $line | cut -d"," -f1)
  
    bam=$(find -L $WD -type f -name "*bam" | fgrep "$sample"| head -1)

    echo $bam

    mkdir -p "$WD"/"$sample"
    cd "$_"
    
    
    ## rather, avoiding running methyldackel multiple times - start     ######

    # bam=$(basename $bam)
    # sample=$(basename $bam .bam)
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

    for interval in $(find $WD/data/ -name "*bed")
    do
        echo $interval
        
        mkdir -p "$WD"/"$sample"/"$(basename $interval .bed)"
        cd "$_"

        mysql --user=genome \
              --host=genome-euro-mysql.soe.ucsc.edu -A -P 3306 \
              -e "select chrom, size from mm9.chromInfo" > mm9.genome

        
        ## if no bedfile with intervals to focus in, proceed
        ##   else, intersect the report with the interval BEDFILE
       

        ## only grepping elements with data (not 0,0 for meth and unmeth reads)
        ## added a column 3 with the CpN position + 1 to be able to run bedtools,
        ##   and took if afterwards with `cut`
        grep -vP '\t0\t0\t' ../"$report" |
            awk '{FS=OFS="\t"; print $1,$2-1,$2,$4,$5,$3,$6,$7}' | \
            "$BEDTOOLS" intersect -a - -b "$interval" -wa > masked
        
        masked_report="$(basename $bam)"_ch_d.cytosine_report_mask_"$(basename $interval)".txt
        mv -f masked "$masked_report"
            

        pigz --processes $nthreads $masked_report
        
        ## mincov
        zcat "$masked_report".gz | awk -v mincov="$mincov" '{
FS=OFS="\t"; 
if ($4+$5 >= mincov)
 print $0
}' | gzip > tmp_covered_"$sample"

        mv -f tmp_covered_"$sample" "$masked_report".gz

        ## motif retrieval
        zcat "$masked_report".gz |
            "$BEDTOOLS" slop -i - \
                        -g mm9.genome \
                        -l 3 -r 4 -s | \
            "$BEDTOOLS" getfasta -fi $MM9 \
                        -bed - \
                        -fo "$masked_report"_cytosine_masked_report_slop.fa \
                        -tab \
                        -s
        
        zcat "$masked_report".gz | paste - "$masked_report"_cytosine_masked_report_slop.fa > tmp_"$sample"

        pigz --processes $nthreads tmp_"$sample"

        rm -rf "$masked_report"_cytosine_masked_report_slop.fa

        ## split into cg and non cg, meth and unmeth in one run
        zcat tmp_"$sample".gz | awk -v cgmeth=cg_meth_"$sample" -v chmeth=ch_meth_"$sample" -v cgunmeth=cg_unmeth_"$sample" -v chunmeth=ch_unmeth_"$sample" '{
OFS=FS="\t"; 
if ($7 == "CG" && $4 > 0) {
  print $4,$5,$6,$7,$8,$9,toupper($10) > cgmeth}
else if ($7 == "CG" && $4 == 0) {
  print $4,$5,$6,$7,$8,$9,toupper($10) > cgunmeth}
else if ($7 != "CG" && $4 > 0) {
  print $4,$5,$6,$7,$8,$9,toupper($10) > chmeth}
else if ($7 != "CG" && $4 == 0) {
  print $4,$5,$6,$7,$8,$9,toupper($10) > chunmeth}
else
  print
}' 

        wc -l cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"    

        # count instances by motif
        for item in  cg_meth_"$sample" cg_unmeth_"$sample" ch_meth_"$sample" ch_unmeth_"$sample"
        do
            cut -f7 "$item" | sort | uniq -c | sed 's/^ *//' > \
                                                   motif_counts_"$item".txt
        done

        rm -f cg_"$sample" ch_"$sample" cg_meth_"$sample" cg_unmeth_"$sample" \
           ch_meth_"$sample" ch_unmeth_"$sample"

        mv -f tmp_"$sample".gz "$sample"_raw_masked_report.txt.gz
        rm -f "$sample"_raw_masked_report.txt.gz
        rm -f "$bam".bai mm9.genome
        
        pigz --processes $nthreads *motif_counts*
        
        ## rather - end                                                     ######
        
        cd $WD/"$sample"
    done            
    
    cd $WD
    
done < "$CONFIG_FILE"

