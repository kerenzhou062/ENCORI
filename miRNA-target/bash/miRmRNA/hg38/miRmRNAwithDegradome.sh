#!/bin/sh

scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/mRNA/align
bedDir=/data/liushun/starbase3/degradome/chow/integration/bed
logDir=$scripts/bash/miRmRNA/log

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/hg38_all_sorted_cons_sb3_align.bed \
  -bed $bedDir/degradome_integration_hg38.bed \
  -output $baseDir/withSeq/hg38_all_sorted_cons_sb3_align.degradome.txt > $logDir/miRmRNAwithDegradome.hg38.log 2>&1
