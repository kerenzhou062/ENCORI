#!/bin/sh

scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/
bedDir=/data/liushun/starbase3/degradome/chow/integration/bed
logDir=$scripts/bash/miRother/log

#circRNA,lncRNA,pseudogene,sncRNA
$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/lncRNA/align/mm10_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_mmu.bed \
  -output $baseDir/lncRNA/align/withSeq/mm10_cons_strict_site_sb3_id.degradome.txt > $logDir/miRlncRNAwithDegradome.mm10.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/circRNA/align/mm10_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_mmu.bed \
  -output $baseDir/circRNA/align/withSeq/mm10_cons_strict_site_sb3_id.degradome.txt > $logDir/miRcircRNAwithDegradome.mm10.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/pseudogene/align/mm10_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_mmu.bed \
  -output $baseDir/pseudogene/align/withSeq/mm10_cons_strict_site_sb3_id.degradome.txt > $logDir/miRpseudogenewithDegradome.mm10.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/sncRNA/align/mm10_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_mmu.bed \
  -output $baseDir/sncRNA/align/withSeq/mm10_cons_strict_site_sb3_id.degradome.txt > $logDir/miRsncRNAwithDegradome.mm10.log 2>&1

