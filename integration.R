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

    stopifnot(names(curr) == names(curr_meth))

    for (motif in names(curr)) {
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

for (item in c('uncovered', 'meth', 'unmeth')) {
    motifs[,item] <- as.numeric(as.character(motifs[,item]))
}


tests <- list()
for (sample in unique(motifs$sample)) {
    tests[[sample]] <- list()
    curr <- motifs[motifs$sample == sample,]
    
    tests[[sample]][['meth_vs_unmeth']] <- chisq.test(curr$meth, curr$unmeth)
    tests[[sample]][['meth_vs_covered']] <- chisq.test(curr$meth, (curr$uncovered + curr$unmeth))
}
