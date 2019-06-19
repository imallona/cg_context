#!/usr/bin/env R
##
## Requires june_2019_wgbs.sh to be run first
##
## 19th june 2019

library(pheatmap)
library(reshape2)
library(ggfortify)
library(heatmaply)

TASK <- "cg_context"
HOME <- '/home/imallona'
TASK <- 'neuro'
WD <-  file.path(HOME, 'cg_context', TASK)
DATA <- file.path(HOME, TASK, 'discretized')

MIN_DEPTH <- 10

samples <-  read.csv(file.path(WD,  'june_2019.conf'), header = FALSE,
                     stringsAsFactors = FALSE)

annot <- data.frame(seq = c('20190524.A-TBWGBS_K',
                            '20190524.A-TBWGBS_B',
                            '20190524.A-TBWGBS_A',
                            '20190524.A-TBWGBS_W'),
                    genotype = c('DNMT3A_KO',
                                 'DNMT3B_addback',
                                 'DNMT3A_addback',
                                 'WT'))

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
for (ssample in gsub('_R1.fastq.gz', '', samples$V1)) {
    curr <- list(cg = data.frame(motif = mdict$cg$motif, row.names = 1),
                 ch = data.frame(motif = mdict$ch$motif, row.names = 1))
    for (context in c('cg', 'ch')) {
        for (thres in c('000', '010', '020', '030', '040', '050', '060', '070',
                        '080', '090', '100')) {
            fn <- file.path(WD, 'discretized',
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
        width = 3000,
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

## fourmers

four <- d

for (ssample in names(d)) {
    for (context in c('cg', 'ch')) {
        ## four[[ssample]][[context]]$fourmer <-substr(, 2, 7)
        ## substr(, 2, 7)
        four[[ssample]][[context]] <- list()
        for (i in 1:ncol(d[[ssample]][[context]])) {
            ## four[[ssample]][[context]][,i] <- tapply(d[[ssample]][[context]][,i],
            ##                                         substr(rownames(d[[ssample]][[context]][i]),
            ##                                                2,7),
            ##                                         function(x) sum(x))

            
            four[[ssample]][[context]][[i]] <- as.data.frame(tapply(d[[ssample]][[context]][,i],
                                                     substr(rownames(d[[ssample]][[context]][i]),
                                                            2,7),
                                                     function(x) sum(x)))
            
            
        }
    }
}

for (ssample in names(four)) {
    for (context in c('cg', 'ch')) {
        foo <- do.call(cbind.data.frame, four[[ssample]][[context]])
        colnames(foo) <- colnames(d[[ssample]][[1]])

        four[[ssample]][[context]] <- foo
        ## four[[context
    }
}


## fourmer processing start

## get proportions (row-wise)

p <- four

for (ssample in names(p)) {
    for (context in c('cg', 'ch')) {
        rowsums <- rowSums(four[[ssample]][[context]])
        for (status in colnames(four[[ssample]][[context]])) {
            four[[ssample]][[context]][,status] <- four[[ssample]][[context]][,status]/rowsums
        }
        rowsums <- NULL
    }
}


## plot start
## let's get a motif sorting for all samples
hc <- list(cg = NULL, ch = NULL)
for (context in c('cg', 'ch')) {
    added <- four[[1]][[context]]
    added[!is.na(added)] <- 0 ## reset
    for (ssample in names(p)) {
        added <- added + four[[ssample]][[context]]
    }
    ## hc[[context]] <- hclust(dist(as.matrix(added[,colSums(added) != 0])), method = 'ward.D2')
    hc[[context]] <- hclust(dist(as.matrix(na.omit(added))), method = 'ward.D2')

}

for (ssample in names(p)) {
    for (context in c('cg', 'ch')) {
        
        png(file.path(WD, sprintf('pheatmap_fourmer_%s_%s.png', context, ssample)),
            width = 1000,
            height = 1800)
        
        pheatmap(as.matrix(four[[ssample]][[context]]),
                 cluster_cols = FALSE,
                 cluster_rows = hc[[context]],
                 main = sprintf('%s %s', context, ssample))
        dev.off()
    }
}



## switch to narrow getting the biggest proportion

m <- four

for (ssample in names(p)) {
    for (context in c('cg', 'ch')) {      
        ## the maximum beta value discrete category for a given motif but
        ## requiring it to greater than 0.05% proportion
        m[[ssample]][[context]]  <- apply(four[[ssample]][[context]],
                     1,
                     function(x) return(tail(as.numeric(names(x[x>0.05])), n = 1)))
    }
}

table(m[[1]][[1]])

## collapsing into dataframes

for (context in c('cg', 'ch')) {
    png(file.path(WD, sprintf('samples_clustering_fourmer_%s.png', context)),
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
    ggsave(gp, filename = file.path(WD, sprintf('samples_pca_fourmer_%s.png', context)),
           width = 5, height = 5, units = "in")
}


## fourmer processing end



## pu vs py plot for each genotype and each methylation status

pu <- c('A', 'G')
py <- c('T', 'C') 

## Taking the four-mers and aggregating


four <- d

for (ssample in names(d)) {
    for (context in c('cg', 'ch')) {
        ## four[[ssample]][[context]]$fourmer <-substr(, 2, 7)
        ## substr(, 2, 7)
        four[[ssample]][[context]] <- list()
        for (i in 1:ncol(d[[ssample]][[context]])) {

            four[[ssample]][[context]][[i]] <- as.data.frame(tapply(d[[ssample]][[context]][,i],
                                                     substr(rownames(d[[ssample]][[context]][i]),
                                                            2,7),
                                                     function(x) sum(x)))
            
            
        }
    }
}

for (ssample in names(four)) {
    for (context in c('cg', 'ch')) {
        foo <- do.call(cbind.data.frame, four[[ssample]][[context]])
        colnames(foo) <- colnames(d[[ssample]][[1]])

        four[[ssample]][[context]] <- foo
  
    }
}

## to aggregate we add a pu or py column to each register

for (ssample in names(four)) {
    four[[ssample]][['cg']][,'chem'] <- NA
    
    four[[ssample]][['cg']][,'chem'] <-  ifelse(
        test = substr(rownames(four[[ssample]][['cg']]),5,5) %in% pu,
        yes = 'pu',
        no = 'py')
    
    stopifnot(all(!is.na(four[[ssample]][['cg']][,'chem'])))
    
}

## aggregating by being pu or py

pupy <- list()
for (ssample in names(four)) {
    pupy[[ssample]] <- aggregate(. ~ chem, four[[ssample]][['cg']], sum)
    rownames(pupy[[ssample]]) <- pupy[[ssample]]$chem
    pupy[[ssample]] <- pupy[[ssample]][,-1]    
}


pupyp <- pupy
for (ssample in names(pupyp)) {

    rowsums <- rowSums(pupyp[[ssample]])
    for (status in colnames(pupyp[[ssample]])) {

        pupyp[[ssample]][,status] <- pupyp[[ssample]][,status]/rowsums
    }
    rowsums <- NULL
    pupyp[[ssample]]$chem <- rownames(pupyp[[ssample]])

}

melted <- melt(pupyp, vars = 'chem')
colnames(melted) <- c('base', 'meth_status', 'proportion_of_motifs', 'genotype')

melted$meth_status <- as.character(melted$meth_status)

melted$meth_status[melted$meth_status == '000'] <- 'unmeth'

melted$meth_status[melted$meth_status == '010'] <- '<10%'
melted$meth_status[melted$meth_status == '020'] <- '10%<x<20%'
melted$meth_status[melted$meth_status == '030'] <- '20%<x<30%'
melted$meth_status[melted$meth_status == '040'] <- '30%<x<40%'
melted$meth_status[melted$meth_status == '050'] <- '40%<x<50%'
melted$meth_status[melted$meth_status == '060'] <- '50%<x<60%'
melted$meth_status[melted$meth_status == '070'] <- '60%<x<70%'
melted$meth_status[melted$meth_status == '080'] <- '70%<x<80%'
melted$meth_status[melted$meth_status == '090'] <- '80%<x<90%'
melted$meth_status[melted$meth_status == '100'] <- '90%<x<100%'
melted$meth_status <- factor(melted$meth_status,
                             levels = c('unmeth',
                                        '<10%',   '10%<x<20%',   '20%<x<30%',
                                        '30%<x<40%',  '40%<x<50%', '50%<x<60%',
                                        '60%<x<70%', '70%<x<80%','80%<x<90%',  '90%<x<100%'))
                                        
mg <- ggplot(melted, aes(x = meth_status,
                         y = proportion_of_motifs,
                         colour = factor(genotype))) + geom_point() +
    xlab('methylation (M/[M+U] reads)') +
    ylab('proportion of motifs')

mg <- mg + facet_grid(. ~ base, margins = TRUE)


ggsave(mg +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) ,
       filename = file.path(WD, sprintf('proportions_cg_1.png')),
        width = 10, height = 5, units = "in")

mg <- ggplot(melted, aes(x = meth_status,
                         y = proportion_of_motifs,
                         colour = factor(base))) + geom_point()  +
    xlab('methylation (M/[M+UM] reads)') +
    ylab('proportion of motifs')

mg <- mg + facet_wrap(. ~ genotype, ncol =2)

ggsave(mg +  theme(axis.text.x = element_text(angle = 90, hjust = 1)),
       filename = file.path(WD, sprintf('proportions_cg_2.png')),
        width = 8, height = 5, units = "in")




## what about plotting proportions according to whether purine    ####
## or pyrimidine, but getting the real motifs?                    ###



for (ssample in names(four)) {
    four[[ssample]][['cg']][,'chem'] <- NA
    
    four[[ssample]][['cg']][,'chem'] <-  ifelse(
        test = substr(rownames(four[[ssample]][['cg']]),5,5) %in% pu,
        yes = 'pu',
        no = 'py')
    
    stopifnot(all(!is.na(four[[ssample]][['cg']][,'chem'])))
    
}


for (ssample in names(four)) {
    four[[ssample]][['cg']][,'chem'] <- NA
    
    four[[ssample]][['cg']][,'chem'] <-  ifelse(
        test = substr(rownames(four[[ssample]][['cg']]),5,5) %in% pu,
        yes = 'pu',
        no = 'py')
    
    stopifnot(all(!is.na(four[[ssample]][['cg']][,'chem'])))
    
}



fourp <- four
for (ssample in names(fourp)) {

    curr <- fourp[[ssample]]$cg
    rowsums <- rowSums(curr[,setdiff(colnames(curr), 'chem')])
    for (status in setdiff(colnames(curr), 'chem')) {

        curr[,status] <- curr[,status]/rowsums
    }
    rowsums <- NULL
    ##curr$chem <- rownames(curr)
    fourp[[ssample]] <- curr
    fourp[[ssample]]$genotype <- ssample
    fourp[[ssample]]$motif <- rownames(curr)

}

fourp <- do.call(rbind.data.frame, fourp)

melted <- melt(fourp, vars = c('chem', 'genotype', 'motif'))

colnames(melted) <- c('base', 'genotype', 'motif',  'meth_status', 'proportion_of_motifs')

melted$meth_status <- as.character(melted$meth_status)

melted$meth_status[melted$meth_status == '000'] <- 'unmeth'

melted$meth_status[melted$meth_status == '010'] <- '<10%'
melted$meth_status[melted$meth_status == '020'] <- '10%<x<20%'
melted$meth_status[melted$meth_status == '030'] <- '20%<x<30%'
melted$meth_status[melted$meth_status == '040'] <- '30%<x<40%'
melted$meth_status[melted$meth_status == '050'] <- '40%<x<50%'
melted$meth_status[melted$meth_status == '060'] <- '50%<x<60%'
melted$meth_status[melted$meth_status == '070'] <- '60%<x<70%'
melted$meth_status[melted$meth_status == '080'] <- '70%<x<80%'
melted$meth_status[melted$meth_status == '090'] <- '80%<x<90%'
melted$meth_status[melted$meth_status == '100'] <- '90%<x<100%'
melted$meth_status <- factor(melted$meth_status,
                             levels = c('unmeth',
                                        '<10%',   '10%<x<20%',   '20%<x<30%',
                                        '30%<x<40%',  '40%<x<50%', '50%<x<60%',
                                        '60%<x<70%', '70%<x<80%','80%<x<90%',  '90%<x<100%'))
                                        


                                     
mg <- ggplot(melted, aes(x = meth_status,
                         y = proportion_of_motifs,
                         colour = factor(genotype))) + geom_boxplot() +
    xlab('methylation (M/[M+U] reads)') +
    ylab('proportion of motifs')

mg <- mg + facet_grid(. ~ base, margins = TRUE)


ggsave(mg +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) ,
       filename = file.path(WD, sprintf('proportions_cg_3.png')),
        width = 15, height = 5, units = "in")

mg <- ggplot(melted, aes(x = meth_status,
                         y = proportion_of_motifs,
                         colour = factor(base))) + geom_boxplot()  +
    xlab('methylation (M/[M+UM] reads)') +
    ylab('proportion of motifs')

mg <- mg + facet_wrap(. ~ genotype, ncol =2)

ggsave(mg +  theme(axis.text.x = element_text(angle = 90, hjust = 1)),
       filename = file.path(WD, sprintf('proportions_cg_4.png')),
        width = 8, height = 5, units = "in")



## getting proportions


## aggregate(four[[ssample]][['cg']]$`000`, by=list(chem = four[[ssample]][['cg']]$chem), FUN=sum)

stop('till here')

## still this maybe should be better normalized by overall dnameth level, or maybe comparing the statuses, more than 0.1 meth, more than 0.2 meth etc? for unmeth and meth statuses


## ## collapse to 4-mers
## cgfour <- cg

## cgfour$sixmer <- rownames(cgfour)
## cgfour$fourmer <- substr(cgfour$sixmer, 2, 7)


## ## for (ssample in samples$V1) {
## ##     foo <- rowsum(cgfour[,1], cgfour$fourmer, reorder = TRUE)
## ## }

## unique(substr(mdict$cg$motif, 2, 7))

## for (ssample in samples$V1) {
##     foo <- as.data.frame(tapply(cgfour[,ssample], cgfour$fourmer, function(x) sum(x)))

##     foo <- rowsum(cgfour[,1], cgfour$fourmer, reorder = TRUE)
## }


## remember to normalize!
