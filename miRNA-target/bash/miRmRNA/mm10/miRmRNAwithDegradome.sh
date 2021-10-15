#!/bin/sh

scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/mRNA/align
bedDir=/data/liushun/starbase3/degradome/chow/integration/bed
logDir=$scripts/bash/miRmRNA/log

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/mm10_all_sorted_cons_sb3_align.bed \
  -bed $bedDir/degradome_integration_mmu.bed \
  -output $baseDir/withSeq/mm10_all_sorted_cons_sb3_align.degradome.txt > $logDir/miRmRNAwithDegradome.mm10.log 2>&1
