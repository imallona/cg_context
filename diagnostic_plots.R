#!/usr/bin/env R
##                                        
## postprocesses bowtie and bwameth mappings (taupo)
##
## 17th apr 2018
## Izaskun Mallona
## GPL 

TASK <- "cg_context_bulk"
HOME <- '/home/imallona'
WD <-  file.path(HOME, TASK)
DATA <- file.path(HOME, 'data')

MIN_DEPTH <- 10

fns <- c('20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked_stranded.txt',
         '20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40_stranded.txt')

fd <- list()
for (fn in fns) {
    toy <- read.table(pipe(sprintf('fgrep -w chr17 %s',
                                        file.path(WD, fn))), header = FALSE)[,c(8,9,4,5,17,18,13,14)]
    ## watson and crick
    colnames(toy) <- c('loc_w', 'seq_w', 'meth_w', 'unmeth_w',
                       'loc_c', 'seq_c', 'meth_c', 'unmeth_c')
    str(toy)

    ## toy <- toy[((toy$meth_w + toy$unmeth_w) >= MIN_DEPTH & (toy$meth_c + toy$unmeth_c) >= MIN_DEPTH),]
    toy$beta_w <- toy$meth_w/(toy$meth_w + toy$unmeth_w)
    toy$beta_c <- toy$meth_c/(toy$meth_c + toy$unmeth_c)

    ## summary(toy)

    ## cor(toy$beta_w, toy$beta_c, method = 'spearman')
    ## png()
    ## plot(toy$beta_w, toy$beta_c)
    ## dev.off()

    fd[[fn]] <- toy
}

names(fd)
names(fd) <- c('3a1_bwa', '3a1_bt2')

table(as.character(fd[['3a1_bwa']]$loc_w) == as.character(fd[['3a1_bt2']]$loc_w))

cor(fd[['3a1_bwa']]$beta_w, fd[['3a1_bt2']]$beta_w, method = 'spearman', use = 'pairwise.complete.obs')
cor(fd[['3a1_bwa']]$beta_c, fd[['3a1_bt2']]$beta_c, method = 'spearman', use = 'pairwise.complete.obs')

set.seed(1)
idx <- sample(x = 1:nrow(fd[[1]]), size = 10000, replace = FALSE)
png(file.path(WD, 'preliminar_mappers.png'), width = 1000, height = 1000 )
par(mfrow = c(2,2))
plot(fd[['3a1_bwa']]$beta_w[idx], fd[['3a1_bt2']]$beta_w[idx], pch = 19)
plot(fd[['3a1_bwa']]$beta_c[idx], fd[['3a1_bt2']]$beta_c[idx], pch = 19)

plot(fd[['3a1_bwa']]$beta_w[idx], fd[['3a1_bwa']]$beta_c[idx], pch = 19)
plot(fd[['3a1_bt2']]$beta_w[idx], fd[['3a1_bt2']]$beta_c[idx], pch = 19)

dev.off()

png(file.path(WD, 'preliminar_mappers_2.png'), width = 1000, height = 1000 )
par(mfrow = c(2,2))
plot(log10(fd[['3a1_bwa']]$meth_w[idx] + 1), log10(fd[['3a1_bt2']]$meth_w[idx] + 1), pch = 19)
plot(log10(fd[['3a1_bwa']]$unmeth_w[idx] + 1), log10(fd[['3a1_bt2']]$unmeth_w[idx] + 1), pch = 19)

plot(log10(fd[['3a1_bwa']]$meth_c[idx] + 1), log10(fd[['3a1_bt2']]$meth_c[idx] + 1), pch = 19)
plot(log10(fd[['3a1_bwa']]$unmeth_c[idx] + 1), log10(fd[['3a1_bt2']]$unmeth_c[idx] + 1), pch = 19)

dev.off()
## not high!
