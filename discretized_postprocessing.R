#!/usr/bin/env R
##
##
## 18th Dec 2018

TASK <- "cg_context"
HOME <- '/home/imallona/mnt/nfs'
WD <-  file.path(HOME, TASK)
DATA <- file.path(HOME, TASK, 'dec_2018')

MIN_DEPTH <- 10

samples <-  read.csv(file.path(WD, 'dec_2018', 'dec_2018.conf'), header = FALSE)

## dictionaries of the whole set of motifs
nucl <- c('T', 'C', 'G', 'A')

mdict <- list('cg' = NULL, 'ch' = NULL)
mdict$cg <- apply(expand.grid('1' = nucl, '2' = nucl, '3' = nucl,
                              '4' = 'C',  '5' = 'G',
                              '6' = nucl, '7' = nucl, '8' = nucl ),
                  1,
                  function(x) paste(x, collapse = ''))



mdict$ch <- apply(expand.grid('1' = nucl, '2' = nucl, '3' = nucl,
                              '4' = 'C',  '5' = c('A', 'T', 'C'),
                              '6' = nucl, '7' = nucl, '8' = nucl),
                  1,
                  function(x) paste(x, collapse = ''))


for (item in names(mdict)) {
    fd <- mdict[[item]]
    fd <- sort(as.character(fd))
    fd <- data.frame(motif = fd, count = 0)
    rownames(fd) <- fd$motif
    mdict[[item]] <- fd
    rm(fd)
}



d <- list()
for (ssample in samples$V1) {
    curr <- list(cg = data.frame(motif = mdict$cg$motif, row.names = 1),
                 ch = data.frame(motif = mdict$ch$motif, row.names = 1))
    for (context in c('cg', 'ch')) {
        for (thres in c('02', '04', '06', '08', '1')) {
            fn <- file.path(WD, 'dec_2018', 'discretized',
                                      sprintf('%s_bwameth_default_motif_counts_%s_%s.txt',
                                              ssample, context, thres))
            if (file.exists(fn)) {
                fd <- read.table(fn, col.names = c('count', 'motif'),
                                  stringsAsFactors = FALSE)

                ## removing Ns and short motifs
                rownames(fd) <- fd$motif
                fd <- fd[grep('N', fd$motif, invert = TRUE),]
                fd <- fd[nchar(fd$motif) == 8,]
                fd <- fd[sort(fd$motif),]

                ## filling with zeroes the motifs that are not present in some dataset

                mis <- setdiff(rownames(mdict[[context]]), rownames(fd))
                fd <- rbind(fd, mdict[[context]][mis,])
                ## lexicographical sorting
                fd <- fd[order(as.character(fd$motif)),]
               
               
                
            } else {
                ## probably is missing due to the fact no motif is fully meth
                fd <- mdict[[context]]
            }

            stopifnot(nrow(fd) == nrow(mdict[[context]]))
            curr[[context]][,thres] <- fd$count
        }
    }
    d[[ssample]] <- curr
}


## get proportions
p <- d

for (ssample in names(p)) {
    for (context in c('cg', 'ch')) {
        rowsums <- rowSums(p[[ssample]][[context]])
        for (status in colnames(p[[ssample]][[context]])) {
            p[[ssample]][[context]][,status] <- p[[ssample]][[context]][,status]/rowsums
        }
        rowsums <- NULL
    }
}

## switch to narrow getting the biggest proportion

m <- p

for (ssample in names(p)) {
    for (context in c('cg', 'ch')) {
        ## m[[ssample]][[context]] <- apply(p[[ssample]][[context]], 1, which.max)
        ## m[[ssample]][[context]] <- apply(p[[ssample]][[context]], 1, function(x) max.col(x, ties.method="first"))
        m[[ssample]][[context]] <- apply(p[[ssample]][[context]],
                                         1,
                                         function(x) names(which.max(x)))
    }
}

table(m[[1]][[1]])

## collapsing into dataframes

cg <- as.data.frame( sapply(m, function(x) return(x$cg)))
ch <- as.data.frame( sapply(m, function(x) return(x$ch)))

## plot and cluster

## still this maybe should be better normalized by overall dnameth level

## remember to normalize!
