#!/bin/bash
##
##
##

TASK="cg_context"
HOME="/home/imallona"
WD="$HOME"/"$TASK"
NFS="$HOME"/mnt/nfs/cg_context

cd $NFS

echo Started at $(date)

/usr/bin/R -e "rmarkdown::render('integration.Rmd')"

echo Ended at $(date)
