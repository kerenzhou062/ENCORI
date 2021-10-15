#!/bin/sh

scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/
bedDir=/data/liushun/starbase3/degradome/chow/integration/bed
logDir=$scripts/bash/miRother/log

#circRNA,lncRNA,pseudogene,sncRNA
$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/lncRNA/align/hg38_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_hg38.bed \
  -output $baseDir/lncRNA/align/withSeq/hg38_cons_strict_site_sb3_id.degradome.txt > $logDir/miRlncRNAwithDegradome.hg38.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/circRNA/align/hg38_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_hg38.bed \
  -output $baseDir/circRNA/align/withSeq/hg38_cons_strict_site_sb3_id.degradome.txt > $logDir/miRcircRNAwithDegradome.hg38.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/pseudogene/align/hg38_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_hg38.bed \
  -output $baseDir/pseudogene/align/withSeq/hg38_cons_strict_site_sb3_id.degradome.txt > $logDir/miRpseudogenewithDegradome.hg38.log 2>&1

$scripts/miRmRNAwithDegradome.py -cpu 30 -input $baseDir/sncRNA/align/hg38_cons_strict_site_sb3_id.bed \
  -bed $bedDir/degradome_integration_hg38.bed \
  -output $baseDir/sncRNA/align/withSeq/hg38_cons_strict_site_sb3_id.degradome.txt > $logDir/miRsncRNAwithDegradome.hg38.log 2>&1

