#!/bin/sh
base="/data/liushun/starbase3"
scripts=$base/scripts
baseDir=$base/miRpredict/mRNA/align
bedDir=$base/clip/rbsSeekerBed
logDir=$scripts/bash/miRmRNA/log

cd $baseDir

$scripts/miRmRNAwithCLIP.py -cpu 30 -input hg38_all_sorted_cons_sb3_align.bed \
  -agoRef $bedDir/reference/hg38_ago.txt \
  -bed $bedDir/formattedBed/hg38 \
  -meta $bedDir/hg38_all_dataset_info_table.txt \
  -output withSeq/hg38_all_sorted_cons_sb3_align.CLIP.txt > $logDir/miRmRNAwithCLIP.hg38.log 2>&1

