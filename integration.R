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
library(qwraps2)

TASK <- "cg_context_bulk"
HOME <- '/home/imallona'
WD <-  file.path(HOME, TASK)
DATA <- file.path(HOME, 'data')

MIN_DEPTH <- 10
NFS <- file.path(HOME, 'mnt', 'nfs')


beta2m <- function(beta) {
    m <- log2(beta/(1 - beta))
    return(m)
}

m2beta <- function(m) {
    beta <- 2^m/(2^m + 1)
    return(beta)
}

setwd(WD)

fns <- file.path(NFS, list.files(NFS, "*stranded.txt.gz", recursive = TRUE))


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

samples_annot <- read.table(text ='sample,genotype,seq
20151223.B-MmES_TKOD3A1c1-3_R,tko+d3a1,bwa_hiseq2k
20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40_stranded.txt,tko+d3a1,bt2_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_6_ACAGTGA_L006_001,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_001,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_002,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16209_131212_SN792_0303_AD2AJ9ACXX_lane6_Undetermined_L006_003,tko+d3a2,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_7_CAGATCA_L007_001,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_001,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_002,tko+d3b1,bwa_hiseq2k
BSSE_QGF_16210_131212_SN792_0303_AD2AJ9ACXX_lane7_Undetermined_L007_003,tko+d3b1,bwa_hiseq2k
SRR2878520,oocyte,bwa_hiseq1k
SRR1274742,tko+d3a2,bwa_hiseq2k
SRR1274743,tko+d3a2,bwa_miseq
SRR1274744,tko+d3b1,bwa_miseq
SRR1274745,tko+d3b1,bwa_hiseq2k
SRR1653162,qko+d3b1,bwa_miseq
SRR2878513,oocyte,bwa_hiseq1k
SRR299053,stadler_es,bwa_hiseq_single
SRR299054,stadler_es,bwa_hiseq_single
SRR299055,stadler_es,bwa_hiseq_single
SRR299056,stadler_es,bwa_gaiix_single
SRR299057,stadler_es,bwa_gaiix_single
SRR299058,stadler_es,bwa_gaiix_single
SRR299059,stadler_es,bwa_gaiix_single
SRR299060,stadler_es,bwa_gaiix_single
SRR299061,stadler_es,bwa_gaiix_single
SRR299062,stadler_es,bwa_gaiix_single', header = TRUE, sep = ',')


rownames(samples_annot) <- samples_annot$sample

stopifnot(length(fd) == nrow(samples_annot))
names(fd)

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


## till here

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
smotifs$sample <- gsub('_bwameth_default_dup_marked_stranded.txt.gz', '', smotifs$sample)

save(smotifs, file = sprintf('stranded_smotifs_%s.RData', format(Sys.time(), "%d_%b_%Y")))


## tabular representation

## normalize by overall enrichment



## norm_meth <- with(smotifs, tapply(meth,
##                             list("sample" = sample, "motif" = short),
##                             function(x) sum(x, na.rm = TRUE)))


## norm_unmeth <- with(smotifs, tapply(unmeth,
##                                     list("sample" = sample, "motif" = short),
##                                     function(x) sum(x, na.rm = TRUE)))


## norm <- norm_meth/(norm_meth + norm_unmeth)

## traspose this

norm_meth <- with(smotifs, tapply(meth,
                            list("motif" = short, "sample" = sample),
                            function(x) sum(x, na.rm = TRUE)))


norm_unmeth <- with(smotifs, tapply(unmeth,
                            list("motif" = short, "sample" = sample),
                            function(x) sum(x, na.rm = TRUE)))

norm <- norm_meth/(norm_meth + norm_unmeth)
colSums(norm, na.rm = TRUE)
## this still is not normalized by sample
## second round normalization, by meth/unmeth marginals

tot_norm_meth <- with(smotifs, tapply(meth,
                            list("sample" = sample),
                            function(x) sum(x, na.rm = TRUE)))

## norm2_meth <- norm_meth/as.vector(tot_norm_meth)


tot_norm_unmeth <- with(smotifs, tapply(unmeth,
                            list("sample" = sample),
                            function(x) sum(x, na.rm = TRUE)))

norm_sample <- tot_norm_meth/(tot_norm_meth + tot_norm_unmeth)


## norm2_unmeth <- norm_unmeth/as.vector(tot_norm_unmeth)

## norm2 <- norm2_meth / (norm2_meth + norm2_unmeth)

normalized <- apply(norm, 1, function(x) return(x/norm_sample))

write.csv(normalized, file = 'normalized.csv')


## plot start
## transform from wide to long for lattice 
normalized_melted <- melt(normalized)


png('normalized_%03d.png', width = 1500, height = 1500)
for (annot in colnames(samples_annot)) {

    print(xyplot(value ~ as.factor(motif) | as.factor(sample),
           data = normalized_melted,
           auto.key = list(columns = 1),
           jitter.x=TRUE,       
           group = samples_annot[normalized_melted$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))


}

dev.off()

## png('normalized_variant_%03d.png', width = 2000, height = 2000)
## for (annot in colnames(samples_annot)) {


##     print(xyplot(value ~ as.factor(sample) | as.factor(motif),
##            data = normalized_melted,
##            auto.key = list(columns = 1, corner = c(1,0)),
##            par.settings = list(layout.widths = list(right.padding = 25)),
##            jitter.x=TRUE,       
##            group = samples_annot[normalized_melted$sample, annot],
##            pch = 19,
##            cex = 0.5,
##            scales=list(x=list(rot=90)),
##            layout = c(6,6)))


## }

## dev.off()




png('normalized_variant_%03d.png', width = 1500, height = 2000)
for (annot in colnames(samples_annot)) {


    print(xyplot(value ~ as.factor(sample) | as.factor(motif),
           data = normalized_melted,
           auto.key = list(columns = 1, corner = c(0.3,1)),
           par.settings = list(layout.widths = list(right.padding = 10)),
           jitter.x=TRUE,       
           group = samples_annot[normalized_melted$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))


}

dev.off()


## plot end


## colSums(norm/norm2, na.rm = TRUE)


## colSums(norm*norm2, na.rm = TRUE)

## this is not what intended, just extract those items that are methylated and compare to the ones represented

'
for sample in fd
   filter those 
'


## plotting

## getting rid of ns
smotifs <- smotifs[grep('n', smotifs$short, invert = TRUE),]

png('new_with_strand_%03d.png', width = 1500, height = 1500)
for (annot in colnames(samples_annot)) {

    print(xyplot(meth_vs_represented_ratio ~ as.factor(paste(short, strand)) | as.factor(sample),
           data = smotifs,
           auto.key = list(columns = 1),
           jitter.y=TRUE,       
           group = samples_annot[smotifs$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(5,5)))

    ## xyplot(meth_vs_represented_ratio ~ as.factor(strand) | as.factor(sample),
    ##        data = smotifs,
    ##        auto.key = list(columns = 1),
    ##        jitter.y=TRUE,       
    ##        group = samples_annot[smotifs$sample, annot],
    ##        pch = 19,
    ##        cex = 0.5,
    ##        scales=list(x=list(rot=90)),
    ##        layout = c(5,5))
}

dev.off()





png('new_with_strand_next_no_jitter_%03d.png', width = 1000, height = 2000)
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

png('new_with_strand_next_jitter_%03d.png', width = 1000, height = 2000)
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


## anyway, needed to be normalized by the methylation level




## some representations on beta values

fns <- file.path(gz_path, list.files(NFS, "*stranded.txt.gz"))

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
betas_df$sample <- gsub('_bwameth_default_dup_marked_stranded.txt.gz', '', betas_df$sample)


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


## transform to mvalue and normalize to the average betavalue


## normalizer <- apply(means[,1:14],2, function(x) sum(na.omit(x))/length(na.omit(x)))


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

dev.off()

## let's print the m values with dots as usuall


betas_strand <- betas_df
betas_strand <- betas_strand[,c('loc_w', 'seq_w', 'meth_w', 'unmeth_w', 'beta_w', 'sample')]
betas_strand$strand <- 'watson'

tmp <- data.frame(betas_df[,c('loc_c', 'seq_c', 'meth_c', 'unmeth_c', 'beta_c', 'sample')],
                  strand = 'crick')

colnames(tmp) <- colnames(betas_strand) <- c('loc', 'seq', 'meth', 'unmeth', 'beta', 'sample', 'strand')
betas_strand <- rbind(betas_strand, tmp)
table(betas_strand$strand)

betas_strand$short <- substr(betas_strand$seq, 2, 5)


betas_df$sample <- gsub('.gz', '', betas_df$sample)
table(betas_df$sample %in% rownames(samples_annot))
unique(betas_df$sample)[!unique(betas_df$sample) %in% rownames(samples_annot)]

betas_strand$sample <- gsub('.gz', '', betas_strand$sample)

png('jitter_m_values_%03d.png', width = 1500, height = 1500)
for (annot in colnames(samples_annot)) {

    print(xyplot(beta2m(beta) ~ as.factor(paste(short, strand)) | as.factor(sample),
           data = betas_strand,
           auto.key = list(columns = 1),
           jitter.x=TRUE,       
           group = samples_annot[betas_strand$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))


}

dev.off()

png('other_jitter_m_values_%03d.png', width = 1000, height = 2000)
for (annot in colnames(samples_annot)) {

    print(xyplot(beta2m(beta) ~ as.factor(sample) | as.factor(strand)*as.factor(short),
           data = betas_strand,
           auto.key = list(columns = 1),
           jitter.x=TRUE,       
           group = samples_annot[betas_strand$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))

}

dev.off()


## and now normalize the m value deviation over the avg m value for each sample

## normalizer <- aggregate(x = beta2m(betas_strand$beta),  as.list(betas_strand$sample),  FUN = mean)
## normalizer <- beta2m(dcast(betas_strand, . ~ sample, mean, value.var = 'beta'))
## with(betas_strand, tapply(beta, as.factor(sample), function(x) mean(beta2m(x), na.rm = TRUE)))
## normalizer <- with(betas_strand, tapply(beta, as.factor(sample), function(x) beta2m(mean(x, na.rm = TRUE))))
normalizer <- with(betas_strand, tapply(beta, as.factor(sample), function(x) mean(x, na.rm = TRUE)))

normalizer.min <-with(betas_strand, tapply(beta, as.factor(sample), function(x) min(x, na.rm = TRUE)))
normalizer.max <-with(betas_strand, tapply(beta, as.factor(sample), function(x) max(x, na.rm = TRUE)))



betas_strand$norm_beta <-  NULL
for (sample in unique(betas_strand$sample)) {
    ## betas_strand[betas_strand$sample == sample, 'norm_beta'] <-  (betas_strand[betas_strand$sample == sample, 'beta'] -
    ##                                                                  normalizer[sample])/normalizer[sample]

    betas_strand[betas_strand$sample == sample, 'norm_beta'] <-  scale(betas_strand[betas_strand$sample == sample, 'beta'])

 
}


png('jitter_norm_beta_values_%03d.png', width = 1500, height = 1500)
for (annot in colnames(samples_annot)) {

    print(xyplot(norm_beta ~ as.factor(paste(short, strand)) | as.factor(sample),
           data = betas_strand,
           auto.key = list(columns = 1),
           jitter.x=TRUE,       
           group = samples_annot[betas_strand$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))


}

dev.off()

png('other_jitter_norm_beta_values_%03d.png', width = 1000, height = 2000)
for (annot in colnames(samples_annot)) {

    print(xyplot(norm_beta~ as.factor(sample) | as.factor(strand)*as.factor(short),
           data = betas_strand,
           auto.key = list(columns = 1),
           jitter.x=TRUE,       
           group = samples_annot[betas_strand$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))

}

dev.off()


## what if scaling to mean 0 variance 1?

## mmm, this are already normalized to mean 0 variance 1
mean(betas_strand[betas_strand$sample == sample, 'norm_beta'], na.rm = TRUE)
sd(betas_strand[betas_strand$sample == sample, 'norm_beta'], na.rm = TRUE)




png('boxplot_other_jitter_norm_beta_values_%03d.png', width = 1000, height = 2000)
for (annot in colnames(samples_annot)) {

    print(bwplot(norm_beta~ as.factor(sample) | as.factor(strand)*as.factor(short),
           data = betas_strand,
           auto.key = list(columns = 1),
           jitter.x=TRUE,       
           group = samples_annot[betas_strand$sample, annot],
           pch = 19,
           cex = 0.5,
           scales=list(x=list(rot=90)),
           layout = c(6,6)))

}

dev.off()



## what about the context of meth cpgs?



## sumtable <-
##   list("meth_vs_represented_ratio" =
##        list("min" = ~ min(meth_vs_represented_ratio),
##             "max" = ~ max(meth_vs_represented_ratio),
##             "mean (sd)" = ~ qwraps2::mean_sd(meth_vs_represented_ratio)),
##        "uncovered" =
##        list("min" = ~ min(uncovered),
##             "max" = ~ max(uncovered),
##             "mean (sd)" = ~ qwraps2::mean_sd(uncovered)),
##        "meth" =
##        list("min" = ~ min(meth),
##             "max" = ~ max(meth),
##             "mean (sd)" = ~ qwraps2::mean_sd(meth))
##        )

## ### Overall

## summary_table(dplyr::group_by(smotifs, short), sumtable)


with(smotifs, tapply(meth_vs_represented_ratio,
                     list("sample" = sample, "motif" = short), mean))
