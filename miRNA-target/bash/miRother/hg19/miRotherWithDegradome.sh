#!/bin/sh

scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/
bedDir=/data/liushun/starbase3/degradome/chow/integration/bed
logDir=$scripts/bash/miRother/log

#circRNA,lncRNA,pseudogene,sncRNA
$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/lncRNA/align/hg19_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_hsa.bed \
  -output $baseDir/lncRNA/align/withSeq/hg19_cons_strict_site_sb3_id.degradome.txt > $logDir/miRlncRNAwithDegradome.hg19.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/circRNA/align/hg19_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_hsa.bed \
  -output $baseDir/circRNA/align/withSeq/hg19_cons_strict_site_sb3_id.degradome.txt > $logDir/miRcircRNAwithDegradome.hg19.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/pseudogene/align/hg19_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_hsa.bed \
  -output $baseDir/pseudogene/align/withSeq/hg19_cons_strict_site_sb3_id.degradome.txt > $logDir/miRpseudogenewithDegradome.hg19.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/sncRNA/align/hg19_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_hsa.bed \
  -output $baseDir/sncRNA/align/withSeq/hg19_cons_strict_site_sb3_id.degradome.txt > $logDir/miRsncRNAwithDegradome.hg19.log 2>&1

