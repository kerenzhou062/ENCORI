#!/bin/sh
base="/data/liushun/starbase3"
scripts=$base/scripts
baseDir=$base/miRpredict/mRNA/align
bedDir=$base/clip/rbsSeekerBed
logDir=$scripts/bash/miRmRNA/log

cd $baseDir

$scripts/miRmRNAwithCLIP.py -cpu 30 -input mm10_all_sorted_cons_sb3_align.bed \
  -agoRef $bedDir/reference/mm10_ago.txt \
  -bed $bedDir/formattedBed/mm10 \
  -meta $bedDir/mm10_all_dataset_info_table.txt \
  -output withSeq/mm10_all_sorted_cons_sb3_align.CLIP.txt > $logDir/miRmRNAwithCLIP.mm10.log 2>&1
