#!/usr/bin/env R
##
##
## 18th Dec 2018

library(pheatmap)
library(reshape2)
library(ggfortify)


TASK <- "cg_context"
HOME <- '/home/imallona/mnt/nfs'
WD <-  file.path(HOME, TASK)
DATA <- file.path(HOME, TASK, 'dec_2018')

MIN_DEPTH <- 10

samples <-  read.csv(file.path(WD, 'dec_2018', 'dec_2018.conf'), header = FALSE)

annot <- data.frame(seq = c('20181207.A-WGBS_IL12',
                            '20181207.A-WGBS_IL14',
                            '20181207.A-WGBS_IL18',
                            '20181207.A-WGBS_IL19'),
                    genotype = c('WT_DNMT3A2_in_TKO',
                                 'QC_DNMT3A2_in_TKO',
                                 'WT_DNMT3A2_in_TKO',
                                 'QC_DNMT3A2_in_TKO'))

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
        for (thres in c('000', '010', '020', '030', '040', '050', '060', '070',
                        '080', '090', '100')) {
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


## get proportions (row-wise)

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


## plot start
## let's get a motif sorting for all samples
hc <- list(cg = NULL, ch = NULL)
for (context in c('cg', 'ch')) {
    added <- p[[1]][[context]]
    added[!is.na(added)] <- 0 ## reset
    for (ssample in names(p)) {
        added <- added + p[[ssample]][[context]]
    }
    ## hc[[context]] <- hclust(dist(as.matrix(added[,colSums(added) != 0])), method = 'ward.D2')
    hc[[context]] <- hclust(dist(as.matrix(na.omit(added))), method = 'ward.D2')

}

for (ssample in names(p)) {
    for (context in c('cg', 'ch')) {
        
        png(file.path(WD, sprintf('pheatmap_%s_%s.png', context, ssample)),
            width = 1000,
            height = 1000)
        
        pheatmap(as.matrix(p[[ssample]][[context]]),
                 cluster_cols = FALSE,
                 cluster_rows = hc[[context]],
                 main = sprintf('%s %s', context, ssample))
        dev.off()
    }
}

## plot end

## switch to narrow getting the biggest proportion

m <- p

for (ssample in names(p)) {
    for (context in c('cg', 'ch')) {
        ## m[[ssample]][[context]] <- apply(p[[ssample]][[context]], 1, which.max)
        ## m[[ssample]][[context]] <- apply(p[[ssample]][[context]], 1, function(x) max.col(x, ties.method="first"))

        ## this would get the unmeth always
        
        ## m[[ssample]][[context]] <- apply(p[[ssample]][[context]],
        ##                                  1,
        ##                                  function(x) names(which.max(x)))

        ## the sort of betavalue with the highest value which is represented
        ## m[[ssample]][[context]]  <- apply(p[[ssample]][[context]],
        ##              1,
        ##              function(x) return(tail(as.numeric(names(x[x!=0])), n = 1)))

        ## the maximum beta value discrete category for a given motif but
        ## requiring it to greater than 0.05% proportion
        m[[ssample]][[context]]  <- apply(p[[ssample]][[context]],
                     1,
                     function(x) return(tail(as.numeric(names(x[x>0.05])), n = 1)))
    }
}

table(m[[1]][[1]])

## collapsing into dataframes

for (context in c('cg', 'ch')) {
    png(file.path(WD, sprintf('samples_clustering_%s.png', context)),
        width = 1000,
        height = 1000)

    
    curr <- as.data.frame( sapply(m, function(x) return(x[[context]])))

    ## plot and cluster

    for (i in 1:ncol(curr)) {
        curr[,i] <- as.numeric(as.character(curr[,i]))
    }

    curr <- t(na.omit(curr))
    ## removing nonvariable motifs
    curr <- curr[, - as.numeric(which(apply(curr, 2, var) == 0))]

    pca <- prcomp(curr, center = TRUE, scale = TRUE)
    ## biplot(pca)
    pheatmap(curr, clustering_method = 'ward', cluster_row = TRUE,
             main = sprintf(context))

    dev.off()

    ## png(file.path(WD, sprintf('samples_pca_%s.png', context)),
    ##     width = 1000,
    ##     height = 1000)
    ## what the hell with the replicates!
    
    gp <- autoplot(pca, data = annot, colour = 'genotype', loadings = FALSE, main = context)
    ## dev.off()
    ggsave(gp, filename = file.path(WD, sprintf('samples_pca_%s.png', context)),
           width = 5, height = 5, units = "in")
}

## the upper is probably not ok, rather than most common meth stats, why not evaluating a sort of beta value?

stop('till here')

## still this maybe should be better normalized by overall dnameth level, or maybe comparing the statuses, more than 0.1 meth, more than 0.2 meth etc? for unmeth and meth statuses


## collapse to 4-mers
cgfour <- cg

cgfour$sixmer <- rownames(cgfour)
cgfour$fourmer <- substr(cgfour$sixmer, 2, 7)


## for (ssample in samples$V1) {
##     foo <- rowsum(cgfour[,1], cgfour$fourmer, reorder = TRUE)
## }

unique(substr(mdict$cg$motif, 2, 7))

for (ssample in samples$V1) {
    foo <- as.data.frame(tapply(cgfour[,ssample], cgfour$fourmer, function(x) median(x)))

    foo <- rowsum(cgfour[,1], cgfour$fourmer, reorder = TRUE)
}


## remember to normalize!
