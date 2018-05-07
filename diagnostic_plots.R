#!/usr/bin/env R
##                                        
## postprocesses bowtie and bwameth mappings (taupo)
##
## 17th apr 2018
## Izaskun Mallona
## GPL

library(kmer)
library(PMCMR)

library(reshape2)
library(lattice)

TASK <- "cg_context_bulk"
HOME <- '/home/imallona'
WD <-  file.path(HOME, TASK)
DATA <- file.path(HOME, 'data')

MIN_DEPTH <- 10


setwd(WD)

fns <- c(file.path(WD, c('20151223.B-MmES_TKOD3A1c1-3_R_bwameth_default_dup_marked_stranded.txt',
         '20151223.B-MmES_TKOD3A1c1-3*_cutadapt_sickle_bismark_bt2_*e.deduplicated_mapq40_stranded.txt')))


fd <- list()
for (fn in fns) {
    toy <- read.table(pipe(sprintf('fgrep -w chr17 %s',
                                         fn)), header = FALSE)[,c(8,9,4,5,17,18,13,14)]
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

    toy$beta <-  (toy$meth_c + toy$meth_w)/(toy$meth_c + toy$unmeth_c + toy$meth_w + toy$unmeth_w)

    fd[[fn]] <- toy
}

names(fd)
names(fd) <- c('3a1_bwa', '3a1_bt2')

## table(as.character(fd[['3a1_bwa']]$loc_w) == as.character(fd[['3a1_bt2']]$loc_w))

## cor(fd[['3a1_bwa']]$beta_w, fd[['3a1_bt2']]$beta_w, method = 'spearman', use = 'pairwise.complete.obs')
## cor(fd[['3a1_bwa']]$beta_c, fd[['3a1_bt2']]$beta_c, method = 'spearman', use = 'pairwise.complete.obs')

## set.seed(1)
## idx <- sample(x = 1:nrow(fd[[1]]), size = 10000, replace = FALSE)
## png(file.path(WD, 'preliminar_mappers.png'), width = 1000, height = 1000 )
## par(mfrow = c(2,2))
## plot(fd[['3a1_bwa']]$beta_w[idx], fd[['3a1_bt2']]$beta_w[idx], pch = 19)
## plot(fd[['3a1_bwa']]$beta_c[idx], fd[['3a1_bt2']]$beta_c[idx], pch = 19)

## plot(fd[['3a1_bwa']]$beta_w[idx], fd[['3a1_bwa']]$beta_c[idx], pch = 19)
## plot(fd[['3a1_bt2']]$beta_w[idx], fd[['3a1_bt2']]$beta_c[idx], pch = 19)

## dev.off()

## png(file.path(WD, 'preliminar_mappers_2.png'), width = 1000, height = 1000 )
## par(mfrow = c(2,2))
## plot(log10(fd[['3a1_bwa']]$meth_w[idx] + 1), log10(fd[['3a1_bt2']]$meth_w[idx] + 1), pch = 19)
## plot(log10(fd[['3a1_bwa']]$unmeth_w[idx] + 1), log10(fd[['3a1_bt2']]$unmeth_w[idx] + 1), pch = 19)

## plot(log10(fd[['3a1_bwa']]$meth_c[idx] + 1), log10(fd[['3a1_bt2']]$meth_c[idx] + 1), pch = 19)
## plot(log10(fd[['3a1_bwa']]$unmeth_c[idx] + 1), log10(fd[['3a1_bt2']]$unmeth_c[idx] + 1), pch = 19)

## dev.off()

## ## let's filter by numreads
## d <- list()
## for (thres in seq(from = 0, to = 50, by = 5)) {
##     thres <- as.character(thres)
##     d[[as.character(thres)]] <- list()
##     for (item in names (fd)) {
##         toy <- fd[[item]]
##         toy <- toy[((toy$meth_w + toy$unmeth_w) >= as.numeric(thres) &
##                         (toy$meth_c + toy$unmeth_c) >= as.numeric(thres)),]
##         d[[as.character(thres)]][[item]] <- toy
##         rm(toy)
##     }

    
##     set.seed(1)
##     if (nrow(d[[thres]][[1]]) >= 10000)
##         idx <- sample(x = 1:nrow(d[[thres]][[1]]), size = 10000, replace = FALSE)
##     else
##         idx <- 1:nrow(d[[thres]][[1]])
    
##     png(file.path(WD, sprintf('preliminar_mappers_thres%s.png', thres)), height = 750, width = 1400)
##     par(mfrow = c(2, 4))

##     par(cex.axis = 1.4,
##         cex.lab = 1.4,
##         cex.main = 1.4,
##         cex.sub = 1.4,
##         pty = "s",
##         mar=c(5.1,4.1,4.1,2.1),
##         oma = c(4, 4, 1, 1))
    
##     plot(d[[thres]][['3a1_bwa']]$beta_w[idx], d[[thres]][['3a1_bt2']]$beta_w[idx], pch = 19)

##     title(sprintf('sample %s min depth %s', '3a1' , thres), outer = TRUE, cex.main = 2)
    

##     plot(d[[thres]][['3a1_bwa']]$beta_c[idx], d[[thres]][['3a1_bt2']]$beta_c[idx], pch = 19)

##     plot(d[[thres]][['3a1_bwa']]$beta_w[idx], d[[thres]][['3a1_bwa']]$beta_c[idx], pch = 19)
##     plot(d[[thres]][['3a1_bt2']]$beta_w[idx], d[[thres]][['3a1_bt2']]$beta_c[idx], pch = 19)


##     plot(log10(d[[thres]][['3a1_bwa']]$meth_w[idx] + 1),
##          log10(d[[thres]][['3a1_bt2']]$meth_w[idx] + 1), pch = 19)
##     plot(log10(d[[thres]][['3a1_bwa']]$unmeth_w[idx] + 1),
##          log10(d[[thres]][['3a1_bt2']]$unmeth_w[idx] + 1), pch = 19)

##     plot(log10(d[[thres]][['3a1_bwa']]$meth_c[idx] + 1),
##          log10(d[[thres]][['3a1_bt2']]$meth_c[idx] + 1), pch = 19)
##     plot(log10(d[[thres]][['3a1_bwa']]$unmeth_c[idx] + 1),
##          log10(d[[thres]][['3a1_bt2']]$unmeth_c[idx] + 1), pch = 19)


##     dev.off()
    
## }

## ## density depths, meth and unmeth
## ## plot(density(df[[]]))


## ## sort by coverage rank and represent

## sorted <- fd

## for (item in names(fd)){
##     sorted[[item]] <- rank(fd[[item]]$meth_w + fd[[item]]$meth_c +
##                                fd[[item]]$unmeth_w + fd[[item]]$unmeth_w, ties.method = 'last')
   
## }


## png(file.path(WD, 'rank_comparison.png'))
## set.seed(1)
## idx <- sample(x = 1:length(sorted[[1]]), size = 10000, replace = FALSE)

## plot(sorted[['3a1_bwa']][idx], sorted[['3a1_bt2']][idx])
## dev.off()


## foo <- data.frame(bwa = sorted[['3a1_bwa']],
##                   bt2 = sorted[['3a1_bt2']])

## foo <- foo[order(foo$bwa),]


## png(file.path(WD, 'rank_comparison_sc.png'))
## smoothScatter(foo)
## dev.off()


## ## overall meth stats

## sapply(fd, function(x) table(x$meth_w > 0))
## sapply(fd, function(x) table(x$unmeth_w > 0))
## sapply(fd, function(x) table(x$unmeth_w > 5))
## sapply(fd, function(x) table(x$meth_w > 5))


## ## anyway, let's get the kmers just in case

## kruskal.test(fd[['3a1_bwa']]$beta_c~as.factor(tolower(as.character(fd[['3a1_bwa']]$seq_c))))


## ## monont
## kruskal.test(fd[['3a1_bwa']]$beta_c,
##              as.factor(substr(tolower(as.character(fd[['3a1_bwa']]$seq_c)), 3, 6)))

## ph <- posthoc.kruskal.nemenyi.test(x = fd[['3a1_bwa']]$beta_c,
##                                    g = as.factor(substr(tolower(as.character(fd[['3a1_bwa']]$seq_c)),
##                                        3, 6)),
##                                    method = 'Tukey')
## ## mononuc
## means <- list()

## for (item in names(fd)) {
##     means[[item]] <- tapply(fd[[item]]$beta_c,
##                             as.factor(substr(tolower(as.character(fd[['3a1_bwa']]$seq_c)), 3, 6)),
##                             function(x) mean(x, na.rm = TRUE))
## }

## if (names((means[[1]])) ==  names((means[[1]])))
##     cor.test(na.omit(as.numeric(means[[1]])), na.omit(as.numeric(means[[2]])), method = 'spearman')


## ## dinuc
## means <- list()

## for (item in names(fd)) {
##     means[[item]] <- tapply(fd[[item]]$beta_c,
##                             as.factor(substr(tolower(as.character(fd[['3a1_bwa']]$seq_c)), 2, 7)),
##                             function(x) mean(x, na.rm = TRUE))
## }

## ## if (names((means[[1]])) ==  names((means[[1]])))
## ##     cor.test(na.omit(as.numeric(means[[1]])), na.omit(as.numeric(means[[2]])), method = 'spearman')


## ## is it the same for other KOs?


## ## further <- c(file.path(WD, c(
## ##     'SRR2878513_bwameth_default_stranded.txt',
## ##     'SRR1274742_bwameth_default_stranded.txt', 'SRR1274743_bwameth_default_stranded.txt',
## ##     'SRR1274744_bwameth_default_stranded.txt', 'SRR1274745_bwameth_default_stranded.txt',
## ##     'SRR1653162_bwameth_default_stranded.txt',
## ##     'SRR2878520_bwameth_default_stranded.txt')),
## ##              file.path(HOME, 'cg_context_new_tuncay', list.files(file.path(HOME, 'cg_context_new_tuncay'),
## ##                                                                  "*stranded.txt", recursive = TRUE)))


further <- c(file.path(WD, c(
    'SRR2878513_bwameth_default_stranded.txt',
    'SRR1274742_bwameth_default_stranded.txt', 'SRR1274743_bwameth_default_stranded.txt',
    'SRR1274744_bwameth_default_stranded.txt', 'SRR1274745_bwameth_default_stranded.txt',
    'SRR1653162_bwameth_default_stranded.txt')),
             file.path(HOME, 'cg_context_new_tuncay', list.files(file.path(HOME, 'cg_context_new_tuncay'),
                                                                 "*stranded.txt", recursive = TRUE)))

for (fn in further) {
    print(fn)
    toy <- read.table(pipe(sprintf('fgrep -w chr17 %s',
                                        file.path(fn))), header = FALSE)[,c(8,9,4,5,17,18,13,14)]
    ## watson and crick
    colnames(toy) <- c('loc_w', 'seq_w', 'meth_w', 'unmeth_w',
                       'loc_c', 'seq_c', 'meth_c', 'unmeth_c')

    toy$beta_w <- toy$meth_w/(toy$meth_w + toy$unmeth_w)
    toy$beta_c <- toy$meth_c/(toy$meth_c + toy$unmeth_c)
    toy$beta <-  (toy$meth_c + toy$meth_w)/(toy$meth_c + toy$unmeth_c + toy$meth_w + toy$unmeth_w)

    fd[[basename(fn)]] <- toy
    print('ended')
}

names(fd) <- gsub('_bwameth_default_stranded.txt', '', basename(names(fd)))



means <- list()


for (item in names(fd)) {
    means[[item]] <- tapply(fd[[item]]$beta_c,
                            as.factor(substr(tolower(as.character(fd[[item]]$seq_c)), 3, 6)),
                                       function(x) mean(x, na.rm = TRUE))
    cnames <- names(means[[item]]) ## this does not change
    means[[item]] <- as.numeric(means[[item]])    
}





## so let's stratify each element by genomic compartement and everything

means <- do.call(cbind.data.frame, means)
means$seq <- cnames

setwd('/home/imallona/cg_context_new_tuncay/')
## png()
## plot(means)
## dev.off()

save(means, file = 'means.RData')

## boxplot(means)


## for (i in 1:ncol(means)){
##     means[,i] <- as.vector(means[,i])
## }
    
melted <- melt(means, id.vars = c("seq"))
xyplot(value ~ as.factor(seq) | variable, data = melted)
xyplot(value ~ variable | as.factor(seq), data = melted, autokey = TRUE)

bwplot(value ~ variable | as.factor(seq), data = melted, autokey = TRUE)

# what if normalizing against the avg meth value?

normalizer <- apply(means[,1:14],2, function(x) sum(na.omit(x))/length(na.omit(x)))

## check
mean(means[,1]/normalizer[1], na.rm = TRUE)

normalized_means <- means
for (i in 1:(ncol(means)-1)) {
     normalized_means[,i] <- means[,i]/normalizer[i]
 }


norm_melted <- melt(normalized_means, id.vars = c("seq"))

norm_melted$variable <- as.factor(gsub('_bwameth_default_stranded.txt', '',
                                       as.character(norm_melted$variable)))

png('normalized.png', width = 600, height = 1200)
bwplot(value ~ variable | as.factor(seq), data = norm_melted, autokey = TRUE,
       scales=list(x=list(rot=90, labels=strtrim(levels(norm_melted$variable), 75))))
## bwplot(value ~ variable | as.factor(seq), data = norm_melted, autokey = TRUE,
##        scales=list(x=list(rot=90)))
dev.off()


## the upper normalization makes no sense, but let's use the raw data


set.seed(1)
idx <- sample(1:nrow(fd[[1]]), size = 5000)

## foo <- do.call(function(x) rbind.data.frame(x[idx,]), fd)
foo <- fd[[1]][1,]
foo$sample <-''
foo <- foo[0,]
for (item in names(fd)) {
    curr <- fd[[item]][idx,]
    curr$sample <- item
    foo <- rbind(foo, curr)
}

save(foo, file = 'foo.RData')


## png('test.png')
## densityplot(beta_w ~ as.factor(sample), groups = as.factor(seq_w), data = foo,
##             plot.points = FALSE, auto.key = FALSE)
## dev.off()



## let's try to get some stats
foo <- foo[(foo$meth_w + foo$unmeth_w + foo$meth_c + foo$unmeth_c) > 5,]
## bwplot(beta_w ~ as.factor(sample), groups = as.factor(seq_w), data = foo[idx,],
##        plot.points = FALSE, auto.key = TRUE)

## mononuc
mono.lm <- lm(foo$beta_c~ as.factor(substr(tolower(as.character(foo$seq_c)), 3, 6))+foo$sample)

di.lm <- lm(foo$beta_c~ as.factor(substr(tolower(as.character(foo$seq_c)), 2, 7))+foo$sample)

anova(mono.lm, di.lm)


foo$dinuc_c <- as.factor(substr(tolower(as.character(foo$seq_c)), 6, 7))

    
png('test_%3d.png', width = 1000, height = 1000)
## xyplot(beta_c~ as.factor(substr(tolower(as.character(foo$seq_c)), 2, 7))| sample,
##        data = foo,
##        scales=list(x=list(rot=90)))

bwplot(beta_c~  sample| dinuc_c,
            data = na.omit(foo),
            scales=list(x=list(rot=90)))


dev.off()


## this is rather stupid, let-s focus: count motif appearance and represented and not

motifs <- list()
    ## data.frame(item = NULL,
    ##            motif = NULL,
    ##            meth = 0,
    ##            unmeth = 0)

##@todo speedup!
for (item in names(fd)) {
    for (motif in unique(substr(tolower(as.character(fd[[1]]$seq_c)), 2, 7))) {
        curr <- fd[[item]][substr(tolower(as.character(fd[[item]]$seq_c)), 2, 7) == motif,]

        curr_covered <- curr[((curr$meth_w + curr$unmeth_w) >= 5),]
        curr_uncovered <- nrow(curr) - nrow(curr_covered)
        curr_meth <- nrow(curr_covered[curr_covered$meth_w > 0,])
        curr_unmeth <-  nrow(curr_covered) - curr_meth

        stopifnot(curr_uncovered + curr_unmeth + curr_meth == nrow(curr))
        motifs[[paste(item, motif)]] <- c(item, motif, curr_uncovered, curr_meth, curr_unmeth)
    }
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
