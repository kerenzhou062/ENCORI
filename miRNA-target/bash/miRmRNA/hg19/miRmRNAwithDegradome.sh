#!/bin/sh

scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/mRNA/align
bedDir=/data/liushun/starbase3/degradome/chow/integration/bed
logDir=$scripts/bash/miRmRNA/log

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/hg19_all_sorted_cons_sb3_align.bed \
  -bed $bedDir/degradome_integration_hsa.bed \
  -output $baseDir/withSeq/hg19_all_sorted_cons_sb3_align.degradome.txt > $logDir/miRmRNAwithDegradome.hg19.log 2>&1
