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
library(Hmisc)

TASK <- "cg_context_bulk"
HOME <- '/home/imallona'
WD <-  file.path(HOME, TASK)
DATA <- file.path(HOME, 'data')

MIN_DEPTH <- 10

beta2m <- function(beta) {
    m <- log2(beta/(1 - beta))
    return(m)
}

m2beta <- function(m) {
    beta <- 2^m/(2^m + 1)
    return(beta)
}




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
                                                             "*stranded.txt.gz", recursive = TRUE)),
         file.path(HOME, 'cg_context', list.files(file.path(HOME, 'cg_context'),
                                                  "*stranded.txt.gz", recursive = TRUE)))


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

## getting the motifs for both strands, and to integrate them

smotifs <- list()
for (item in names(fd)) {
    ## for (motif in unique(substr(tolower(as.character(fd[[1]]$seq_c)), 2, 7))) {


    
    fd[[item]]$factor_w <- substr(tolower(as.character(fd[[item]]$seq_w)), 2, 7)
    fd[[item]]$factor_c <- substr(tolower(as.character(fd[[item]]$seq_c)), 2, 7)
    
    curr_covered_w <- with(fd[[item]], tapply(meth_w + unmeth_w, factor_w, function(x) x >= 5))
    curr_meth_w <-  with(fd[[item]], tapply(meth_w , factor_w, function(x) x > 0))

    curr_covered_c <- with(fd[[item]], tapply(meth_c + unmeth_c, factor_c, function(x) x >= 5))
    curr_meth_c <-  with(fd[[item]], tapply(meth_c , factor_c, function(x) x > 0))

    stopifnot(names(curr_covered_w) == names(curr_meth_w))
    stopifnot(names(curr_covered_c) == names(curr_meth_c))

    for (motif in names(curr_covered_w)) {
        ## vector of
        ## sample
        ## motif
        ## cgs uncovered
        ## cgs covered and methylated
        ## cgs covered and unmethylated
        smotifs[[paste(item, motif, 'watson')]] <- c(
            item,
            'watson',
            motif,
            nrow(fd[[item]][fd[[item]]$factor_w == motif,]) -
                sum(curr_covered_w[[motif]]),
            sum(curr_covered_w[[motif]] & curr_meth_w[[motif]]),
            sum(curr_covered_w[[motif]] & ! curr_meth_w[[motif]]))
        
        stopifnot(as.numeric(smotifs[[paste(item, motif)]][3]) +
                      as.numeric(smotifs[[paste(item, motif)]][4]) +
                          as.numeric(smotifs[[paste(item, motif)]][5]) ==
                              nrow(fd[[item]][fd[[item]]$factor_w == motif,]))

    }
    ## this might sound stupid because watson and crick are simetrical and not independent
    ## but anyway, they might end being like that
    for (motif in names(curr_covered_c)) {
        ## vector of
        ## sample
        ## motif
        ## cgs uncovered
        ## cgs covered and methylated
        ## cgs covered and unmethylated
        smotifs[[paste(item, motif, 'crick')]] <- c(
            item,
            'crick',
            motif,
            nrow(fd[[item]][fd[[item]]$factor_c == motif,]) -
                sum(curr_covered_c[[motif]]),
            sum(curr_covered_c[[motif]] & curr_meth_c[[motif]]),
            sum(curr_covered_c[[motif]] & ! curr_meth_c[[motif]]))
        
        stopifnot(as.numeric(smotifs[[paste(item, motif)]][3]) +
                      as.numeric(smotifs[[paste(item, motif)]][4]) +
                          as.numeric(smotifs[[paste(item, motif)]][5]) ==
                              nrow(fd[[item]][fd[[item]]$factor_c == motif,]))

    }
    
                 
    ## }
}


save(smotifs, file = sprintf('stranded_motifs_list_%s.RData', format(Sys.time(), "%d_%b_%Y")))


smotifs <- as.data.frame(do.call(rbind.data.frame, smotifs))
colnames(smotifs) <- c('sample', 'strand', 'motif', 'uncovered', 'meth', 'unmeth')
smotifs$sample <- basename(as.character(smotifs$sample))

for (item in c('uncovered', 'meth', 'unmeth')) {
    smotifs[,item] <- as.numeric(as.character(smotifs[,item]))
}


smotifs$meth_vs_represented_ratio <- smotifs$meth  / (smotifs$meth + smotifs$unmeth)


smotifs$motif <- tolower(as.character(smotifs$motif))

smotifs$short <- substr(smotifs$motif, 2,5)

smotifs$sample <- gsub('_bwameth_default_stranded.txt.gz', '', smotifs$sample)

save(smotifs, file = sprintf('stranded_smotifs_%s.RData', format(Sys.time(), "%d_%b_%Y")))


## plotting

## getting rid of ns
smotifs <- smotifs[grep('n', smotifs$short, invert = TRUE),]

png('with_strand_%03d.png', width = 1500, height = 1500)

xyplot(meth_vs_represented_ratio ~ as.factor(short) | as.factor(sample),
       data = smotifs,
       auto.key = list(columns = 1),
       jitter.y=TRUE,       
       group = samples_annot[smotifs$sample, annot],
       pch = 19,
       cex = 0.5,
       scales=list(x=list(rot=90)),
       layout = c(4,4))

xyplot(meth_vs_represented_ratio ~ as.factor(strand) | as.factor(sample),
       data = smotifs,
       auto.key = list(columns = 1),
       jitter.y=TRUE,       
       group = samples_annot[smotifs$sample, annot],
       pch = 19,
       cex = 0.5,
       scales=list(x=list(rot=90)),
       layout = c(4,4))




dev.off()

png('with_strand_next_no_jitter_%03d.png', width = 1000, height = 2000)
for (annot in colnames(samples_annot)) {
    print(xyplot(meth_vs_represented_ratio ~ as.factor(sample) | as.factor(strand)*as.factor(short),
           data = smotifs,
           auto.key = list(columns = 1),
           jitter.y=TRUE,       
           group = samples_annot[smotifs$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))
}

dev.off()

png('with_strand_next_jitter_%03d.png', width = 1000, height = 2000)
for (annot in colnames(samples_annot)) {
    print(xyplot(meth_vs_represented_ratio ~ as.factor(sample) | as.factor(strand)*as.factor(short),
           data = smotifs,
           auto.key = list(columns = 1),
           jitter.y=FALSE,       
           group = samples_annot[smotifs$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))
}

dev.off()


## some representations on beta values

gz_path <- '/home/imallona/mnt/baubec/imallona2mmanzo/'
fns <- file.path(gz_path, list.files(gz_path, "*stranded.txt.gz"))

MIN_DEPTH <- 5

betas <- list()
for (fn in fns) {
    print(fn)
    toy <- read.table(pipe(sprintf('zcat %s | fgrep -w chr17',
                                   fn)), header = FALSE,
                      stringsAsFactors = FALSE)[,c(8,9,4,5,17,18,13,14)]
    ## watson and crick
    colnames(toy) <- c('loc_w', 'seq_w', 'meth_w', 'unmeth_w',
                       'loc_c', 'seq_c', 'meth_c', 'unmeth_c')
    str(toy)

    

    toy <- toy[((toy$meth_w + toy$unmeth_w) >= MIN_DEPTH | (toy$meth_c + toy$unmeth_c) >= MIN_DEPTH),]
    toy$beta_w <- toy$meth_w/(toy$meth_w + toy$unmeth_w)
    toy$beta_c <- toy$meth_c/(toy$meth_c + toy$unmeth_c)

    toy$beta <-  (toy$meth_c + toy$meth_w)/(toy$meth_c + toy$unmeth_c + toy$meth_w + toy$unmeth_w)

    
    betas[[fn]] <- toy
    betas[[fn]]$sample <- gsub('_bwameth_default_stranded.txt', '', basename(fn))
    print(dim(toy))
}



betas_df <- do.call(rbind.data.frame, betas)

save(betas_df, file = sprintf('stranded_betas_%s.RData', format(Sys.time(), "%d_%b_%Y")))


betas_df$sample <- gsub('_bwameth_default_stranded.txt.gz', '', basename(betas_df$sample))

rownames(betas_df) <- 1:nrow(betas_df)

for (item in c('meth_w', 'unmeth_w', 'meth_c', 'unmeth_c'))  {
    betas_df[,item] <- as.numeric(as.character(betas_df[,item]))
}

for (item in c('seq_w', 'seq_c')) {
    
    betas_df[,item] <- tolower(as.character(betas_df[,item]))

    betas_df[,sprintf('%s_short', item)] <- substr(betas_df[,item], 3,6)
}

betas_df <- betas_df[grep('n', betas_df$seq_w_short, invert = TRUE),]
betas_df <- betas_df[grep('n', betas_df$seq_c_short, invert = TRUE),]

save(betas_df, file = sprintf('stranded_betas_df_%s.RData', format(Sys.time(), "%d_%b_%Y")))



## get some short motifs for the betas

pal <- RColorBrewer::brewer.pal(6, "Set3")

png('violins_betas_%003d.png', width = 1000, height = 2000)

print(bwplot(beta_w ~ as.factor(sample) | as.factor(seq_w_short) ,
       data = betas_df,
       auto.key = list(columns = 4),
       group = samples_annot[betas_df$sample, annot],
             scales=list(x=list(rot=90)),
       panel = function(..., box.ratio) {
           panel.violin(..., col = "transparent",
                        varwidth = FALSE, box.ratio = box.ratio)
           panel.bwplot(..., fill = NULL, box.ratio = .1)
       } ))

print(bwplot(beta_w ~ as.factor(sample) | as.factor(seq_w_short) ,
       data = betas_df,
       auto.key = list(columns = 4),
       group = samples_annot[betas_df$sample, annot],
             scales=list(x=list(rot=90)),
       panel = function(..., box.ratio) {
           panel.violin(..., col = pal,
                        varwidth = FALSE, box.ratio = box.ratio)
           panel.bwplot(..., fill = NULL, box.ratio = .1)
       } ))


print(bwplot(beta_w ~ as.factor(sample) | as.factor(seq_w_short) ,
             data = betas_df,             
             group = samples_annot[betas_df$sample, annot],
             panel = panel.superpose,
             panel.groups = panel.violin),
      col = pal,
      scales=list(x=list(rot=90)))

dev.off()


png('violins_m_%003d.png', width = 1000, height = 2000)

print(bwplot(beta2m(beta_w) ~ as.factor(sample) | as.factor(seq_w_short) ,
             data = betas_df,
             auto.key = list(columns = 4),
             group = samples_annot[betas_df$sample, annot],
             scales=list(x=list(rot=90)),
             panel = function(..., box.ratio) {
                 panel.violin(..., col = "transparent",
                              varwidth = FALSE, box.ratio = box.ratio)
                 panel.bwplot(..., fill = NULL, box.ratio = .1)
             } ))

print(bwplot(beta2m(beta_w) ~ as.factor(sample) | as.factor(seq_w_short) ,
       data = betas_df,
       auto.key = list(columns = 4),
       group = samples_annot[betas_df$sample, annot],
             scales=list(x=list(rot=90)),
       panel = function(..., box.ratio) {
           panel.violin(..., col = pal,
                        varwidth = FALSE, box.ratio = box.ratio)
           panel.bwplot(..., fill = NULL, box.ratio = .1)
       } ))


print(bwplot(beta2m(beta_w) ~ as.factor(sample) | as.factor(seq_w_short) ,
             data = betas_df,             
             group = samples_annot[betas_df$sample, annot],
             panel = panel.superpose,
             panel.groups = panel.violin),
      col = pal,
      scales=list(x=list(rot=90)))


print(bwplot(beta2m(beta_w) ~ as.factor(sample) | as.factor(seq_w_short) ,
             data = betas_df,
             auto.key = list(columns = 4),
             group = samples_annot[betas_df$sample, annot],
             scales=list(x=list(rot=90)),
             panel = function(..., box.ratio) {
                 panel.violin(..., col = "transparent",
                              varwidth = FALSE, box.ratio = box.ratio)
                 panel.bwplot(..., fill = NULL, box.ratio = .1)
                 panel.abline(v=quantile(beta2m(na.omit(betas_df$beta_w)),.5), col.line="red") 
             } ))




## print(densityplot(beta2m(beta_w) ~ as.factor(sample) | as.factor(seq_w_short) ,
##                   data = betas_df,             
##                   group = samples_annot[betas_df$sample, annot],
##                   panel=function(x,...){
##                       panel.densityplot(x,...)
##                       panel.abline(v=quantile(x,.5), col.line="red") 
##                   }))



dev.off()


## ## tests
# does not work, takes huge time to compute
## png('bpplot_m_%003d.png', width = 1000, height = 2000)

## # same as previous but add a spike to give 0.95 interval
## ## bwplot(g ~ x, panel=panel.bpplot, probs=c(.025,seq(.25,.49,by=.01)))
## ## print(bwplot(beta2m(beta_w) ~ as.numeric(sample),
## ##              data = betas_df,
## ##              auto.key = list(columns = 4),
## ##              group = samples_annot[betas_df$sample, annot],
## ##              scales=list(x=list(rot=90)),
## ##              panel=panel.bpplot, probs=c(.025,seq(.25,.49,by=.01))))

## print(bwplot(beta2m(beta_w) ~ as.numeric(sample) | as.factor(seq_w_short),
##              data = betas_df,
##              auto.key = list(columns = 4),
##              group = samples_annot[betas_df$sample, annot],
##              scales=list(x=list(rot=90)),
##              panel=panel.bpplot, probs=c(.025,seq(.25,.49,by=.01))))
## dev.off()
