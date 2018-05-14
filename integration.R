#!/usr/bin/env R
##                                        
## postprocesses bowtie and bwameth mappings (taupo)
##
## 17th apr 2018
## Izaskun Mallona
## GPL

## library(kmer)
## library(PMCMR)

library(reshape2)
library(lattice)

TASK <- "cg_context_bulk"
HOME <- '/home/imallona'
WD <-  file.path(HOME, TASK)
DATA <- file.path(HOME, 'data')

MIN_DEPTH <- 10


setwd(WD)

fns <- c(file.path(WD, c('20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked_stranded.txt.gz',
                         '20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40_stranded.txt.gz')),
         file.path(WD, c(
             ## 'SRR2878520_bwameth_default_stranded.txt.gz',
             'SRR2878513_bwameth_default_stranded.txt.gz',
             'SRR1274742_bwameth_default_stranded.txt.gz',
             'SRR1274743_bwameth_default_stranded.txt.gz',
             'SRR1274744_bwameth_default_stranded.txt.gz',
             'SRR1274745_bwameth_default_stranded.txt.gz',
             'SRR1653162_bwameth_default_stranded.txt.gz')),
         file.path(HOME, 'cg_context_new_tuncay', list.files(file.path(HOME, 'cg_context_new_tuncay'),
                                                             "*stranded.txt.gz", recursive = TRUE)))


## fns_annot <- list(
##     'foo' = 'tko',
    
##     )

fd <- list()
for (fn in fns) {
    toy <- read.table(pipe(sprintf('zcat %s | fgrep -w chr17',
                                         fn)), header = FALSE)[,c(8,9,4,5,17,18,13,14)]
    ## watson and crick
    colnames(toy) <- c('loc_w', 'seq_w', 'meth_w', 'unmeth_w',
                       'loc_c', 'seq_c', 'meth_c', 'unmeth_c')
    str(toy)

    ## toy <- toy[((toy$meth_w + toy$unmeth_w) >= MIN_DEPTH & (toy$meth_c + toy$unmeth_c) >= MIN_DEPTH),]
    toy$beta_w <- toy$meth_w/(toy$meth_w + toy$unmeth_w)
    toy$beta_c <- toy$meth_c/(toy$meth_c + toy$unmeth_c)

    toy$beta <-  (toy$meth_c + toy$meth_w)/(toy$meth_c + toy$unmeth_c + toy$meth_w + toy$unmeth_w)

    fd[[fn]] <- toy
}

## getting the meth tables, just contingency

## motifs <- list()
    ## data.frame(item = NULL,
    ##            motif = NULL,
    ##            meth = 0,
    ##            unmeth = 0)

##@todo speedup!
print('beware this is only watson')

## for (item in names(fd)) {
##     for (motif in unique(substr(tolower(as.character(fd[[1]]$seq_c)), 2, 7))) {
##         curr <- fd[[item]][substr(tolower(as.character(fd[[item]]$seq_c)), 2, 7) == motif,]

##         curr_covered <- curr[((curr$meth_w + curr$unmeth_w) >= 5),]
##         curr_uncovered <- nrow(curr) - nrow(curr_covered)
##         curr_meth <- nrow(curr_covered[curr_covered$meth_w > 0,])
##         curr_unmeth <-  nrow(curr_covered) - curr_meth

##         stopifnot(curr_uncovered + curr_unmeth + curr_meth == nrow(curr))
##         motifs[[paste(item, motif)]] <- c(item, motif, curr_uncovered, curr_meth, curr_unmeth)
##     }
## }



## dat <-data.frame(factor=sample(c("a","b","c"), 10, T), value=rnorm(10))
## r1<-with(dat, tapply(value, factor, mean))
## r1
## r1[["a"]]


motifs <- list()
for (item in names(fd)) {
    ## for (motif in unique(substr(tolower(as.character(fd[[1]]$seq_c)), 2, 7))) {
    fd[[item]]$factor <- substr(tolower(as.character(fd[[item]]$seq_c)), 2, 7)
    curr_covered <- with(fd[[item]], tapply(meth_w + unmeth_w, factor, function(x) x >= 5))

    curr_meth <-  with(fd[[item]], tapply(meth_w , factor, function(x) x > 0))

    stopifnot(names(curr_covered) == names(curr_meth))

    for (motif in names(curr_covered)) {
        ## vector of
        ## sample
        ## motif
        ## cgs uncovered
        ## cgs covered and methylated
        ## cgs covered and unmethylated
        motifs[[paste(item, motif)]] <- c(item,
                                          motif,
                                          nrow(fd[[item]][fd[[item]]$factor == motif,]) -
                                              sum(curr_covered[[motif]]),
                                          sum(curr_covered[[motif]] & curr_meth[[motif]]),
                                          sum(curr_covered[[motif]] & ! curr_meth[[motif]]))

        stopifnot(as.numeric(motifs[[paste(item, motif)]][3]) +
                      as.numeric(motifs[[paste(item, motif)]][4]) +
                          as.numeric(motifs[[paste(item, motif)]][5]) ==
                              nrow(fd[[item]][fd[[item]]$factor == motif,]))

    }
                 
    ## }
}


motifs <- as.data.frame(do.call(rbind.data.frame, motifs))
colnames(motifs) <- c('sample', 'motif', 'uncovered', 'meth', 'unmeth')
motifs$sample <- basename(as.character(motifs$sample))

for (item in c('uncovered', 'meth', 'unmeth')) {
    motifs[,item] <- as.numeric(as.character(motifs[,item]))
}

save(motifs, file = 'motifs.RData')

tests <- list()
for (sample in unique(motifs$sample)) {
    tests[[sample]] <- list()
    curr <- motifs[motifs$sample == sample,]
    
    tests[[sample]][['meth_vs_unmeth']] <- chisq.test(curr$meth, curr$unmeth)
    tests[[sample]][['meth_vs_covered']] <- chisq.test(curr$meth, (curr$uncovered + curr$unmeth))
}

motifs$ratio_m_represented <- motifs$meth  / (motifs$meth + motifs$unmeth)


motifs$motif <- tolower(as.character(motifs$motif))

motifs$short <- substr(motifs$motif, 3,6)

motifs$sample <- gsub('_bwameth_default_stranded.txt.gz', '', motifs$sample)

bwplot(ratio_m_represented ~ sample | as.factor(short),
       data = motifs,
       autokey = TRUE)


bwplot(ratio_m_represented ~ as.factor(short) | as.factor(sample),
       data = motifs,
       autokey = TRUE,
       scales=list(x=list(rot=90)))

dotplot(ratio_m_represented ~ as.factor(short) | as.factor(sample),
       data = motifs,
       autokey = TRUE,
        scales=list(x=list(rot=90)))

## let's try to set up this
## samples_annot <- read.table(text ='sample','genotype','seq'
## '20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked_stranded.txt.gz','tko+d3a1','bwa_hiseq2k'
## '20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40_stranded.txt.gz','tko+d3a1','bt2_hiseq2k'
## 'BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_001','tko+d3a2','bwa_hiseq2k'
## 'BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_001','tko+d3a2','bwa_hiseq2k'
## 'BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_002','tko+d3a2','bwa_hiseq2k'
## 'BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_003','tko+d3a2','bwa_hiseq2k'
## 'BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_001','tko+d3b1','bwa_hiseq2k'
## 'BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_001','tko+d3b1','bwa_hiseq2k'
## 'BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_002','tko+d3b1','bwa_hiseq2k'
## 'BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_003','tko+d3b1','bwa_hiseq2k'
## 'SRR1274742','bwa_hiseq2k','tko+3a2'
## 'SRR1274743','bwa_miseq','tko+3a2'
## 'SRR1274744','bwa_miseq','tko+3b1'
## 'SRR1274745','bwa_hiseq2k','tko+3b1'
## 'SRR1653162','bwa_miseq','qko+3b1'
##                             'SRR2878513','bwa_hiseq1k','oocyte', header = TRUE, sep = ',')

samples_annot <- read.table(text ='sample,genotype,seq
20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked_stranded.txt.gz,tko+d3a1,bwa_hiseq2k
20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40_stranded.txt.gz,tko+d3a1,bt2_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_001,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_001,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_002,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_003,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_001,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_001,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_002,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_003,tko+d3b1,bwa_hiseq2k
SRR1274742,tko+d3a2,bwa_hiseq2k
SRR1274743,tko+d3a2,bwa_miseq
SRR1274744,tko+d3b1,bwa_miseq
SRR1274745,tko+d3b1,bwa_hiseq2k
SRR1653162,qko+d3b1,bwa_miseq
SRR2878513,oocyte,bwa_hiseq1k', header = TRUE, sep = ',')

rownames(samples_annot) <- samples_annot$sample

motifs$color <- c('a', 'b')

## https://www.stat.ubc.ca/~jenny/STAT545A/block16_colorsLatticeQualitative.html

xyplot(ratio_m_represented ~ as.factor(short) | as.factor(sample),
       data = motifs,
       autokey = TRUE,
       jitter.y=TRUE,
       group = motifs$color,
       pch = 19,
       cex = 0.5,
       scales=list(x=list(rot=90)))



for (annot in colnames(samples_annot)) {
    png(sprintf('test_%s.png', annot), width = 1000, height = 2000)

    print(xyplot(ratio_m_represented ~ as.factor(sample) | as.factor(short),
                 data = motifs,
                 auto.key = list(columns = 4),
                 jitter.y=TRUE,
                 group = samples_annot[motifs$sample, annot],
                 pch = 19,
                 cex = 0.5,
                 scales=list(x=list(rot=90)),
                 layout=c(5,4)))

    dev.off()
}

png(sprintf('test_superposed.png'), width = 1000, height = 2000)

print(xyplot(ratio_m_represented ~ as.factor(sample) | as.factor(short),
             data = motifs,
             auto.key = list(columns = 4),
             jitter.y=TRUE,
             group = samples_annot[motifs$sample, 'genotype'],
             pch =  as.numeric(samples_annot[motifs$sample, 'seq']),
             cex = 0.5,
             scales=list(x=list(rot=90)),
             layout=c(5,4)))

dev.off()
