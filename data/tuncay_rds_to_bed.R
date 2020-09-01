#!/usr/bin/env R
##
## Reads Tuncay's GR RDSs placed in folder `./data` with name `*rds` and writes BED files
## Assumed mm9, metadata not found
##
## Izaskun Mallona
## 4th August 2020

path <- 'data'

for (fn in list.files(path, pattern = '.*rds$')) {
    (gr <- readRDS(file.path('data', fn)))

    fd <- data.frame(seqnames = seqnames(gr),
                     starts = start(gr)-1,
                     ends = end(gr),
                     names = c(rep(".", length(gr))),
                     scores = c(rep(".", length(gr))),
                     strands = strand(gr))

    write.table(fd, file = file.path('data', sprintf('%s.bed', gsub('.gr.rds', '', basename(fn)))),
                quote = FALSE, sep = "\t", row.names = F, col.names = FALSE)

    rm(gr)
}


