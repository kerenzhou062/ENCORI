#!/bin/sh
base="/data/liushun/starbase3"
scripts=$base/scripts
baseDir=$base/miRpredict
bedDir=$base/clip/rbsSeekerBed
logDir=$scripts/bash/miRother/log

#circRNA,lncRNA,pseudogene,sncRNA
$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/lncRNA/align/hg19_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/hg19_ago.txt \
  -bed $bedDir/formattedBed/hg19 \
  -meta $bedDir/hg19_all_dataset_info_table.txt \
  -output $baseDir/lncRNA/align/withSeq/hg19_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRlncRNAwithCLIP.hg19.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/circRNA/align/hg19_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/hg19_ago.txt \
  -bed $bedDir/formattedBed/hg19 \
  -meta $bedDir/hg19_all_dataset_info_table.txt \
  -output $baseDir/circRNA/align/withSeq/hg19_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRcircRNAwithCLIP.hg19.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/pseudogene/align/hg19_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/hg19_ago.txt \
  -bed $bedDir/formattedBed/hg19 \
  -meta $bedDir/hg19_all_dataset_info_table.txt \
  -output $baseDir/pseudogene/align/withSeq/hg19_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRpseudogenewithCLIP.hg19.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/sncRNA/align/hg19_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/hg19_ago.txt \
  -bed $bedDir/formattedBed/hg19 \
  -meta $bedDir/hg19_all_dataset_info_table.txt \
  -output $baseDir/sncRNA/align/withSeq/hg19_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRsncRNAwithCLIP.hg19.log 2>&1

