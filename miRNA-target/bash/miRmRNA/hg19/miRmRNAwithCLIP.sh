#!/bin/sh
base="/data/liushun/starbase3"
scripts=$base/scripts
baseDir=$base/miRpredict/mRNA/align
bedDir=$base/clip/rbsSeekerBed
logDir=$scripts/bash/miRmRNA/log

cd $baseDir

$scripts/miRmRNAwithCLIP.py -cpu 30 -input hg19_all_sorted_cons_sb3_align.bed \
  -agoRef $bedDir/reference/hg19_ago.txt \
  -bed $bedDir/formattedBed/hg19 \
  -meta $bedDir/hg19_all_dataset_info_table.txt \
  -output withSeq/hg19_all_sorted_cons_sb3_align.CLIP.txt > $logDir/miRmRNAwithCLIP.hg19.log 2>&1

