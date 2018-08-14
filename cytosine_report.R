#!/usr/bin/env R
##                                        
## postprocesses the cpg meth reports from methyldakel/bismark meth extract
## requires cytosine_report.sh to be run first
##
## 17th apr 2018
## Izaskun Mallona
## GPL 

TASK <- "cytosine_report"
HOME <- '/home/imallona'
WD <-  file.path(HOME, TASK)
DATA <- file.path(HOME, 'data')

MIN_DEPTH <- 10

fn <- '20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked_stranded.txt'
toy <- read.table(pipe(sprintf('fgrep chr10 %s',
                               file.path(WD, fn))), header = FALSE)[,c(8,9,4,5,17,18,13,14)]
## watson and crick
colnames(toy) <- c('loc_w', 'seq_w', 'meth_w', 'unmeth_w',
                   'loc_c', 'seq_c', 'meth_c', 'unmeth_c')
str(toy)

toy <- toy[((toy$meth_w + toy$unmeth_w) >= MIN_DEPTH & (toy$meth_c + toy$unmeth_c) >= MIN_DEPTH),]
toy$beta_w <- toy$meth_w/(toy$meth_w + toy$unmeth_w)
toy$beta_c <- toy$meth_c/(toy$meth_c + toy$unmeth_c)

summary(toy)

cor(toy$beta_w, toy$beta_c, method = 'spearman')
png()
plot(toy$beta_w, toy$beta_c)
dev.off()

## not high!
