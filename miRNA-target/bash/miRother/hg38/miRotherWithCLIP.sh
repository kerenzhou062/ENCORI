#!/bin/sh
base="/data/liushun/starbase3"
scripts=$base/scripts
baseDir=$base/miRpredict
bedDir=$base/clip/rbsSeekerBed
logDir=$scripts/bash/miRother/log

#circRNA,lncRNA,pseudogene,sncRNA
$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/lncRNA/align/hg38_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/hg38_ago.txt \
  -bed $bedDir/formattedBed/hg38 \
  -meta $bedDir/hg38_all_dataset_info_table.txt \
  -output $baseDir/lncRNA/align/withSeq/hg38_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRlncRNAwithCLIP.hg38.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/circRNA/align/hg38_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/hg38_ago.txt \
  -bed $bedDir/formattedBed/hg38 \
  -meta $bedDir/hg38_all_dataset_info_table.txt \
  -output $baseDir/circRNA/align/withSeq/hg38_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRcircRNAwithCLIP.hg38.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/pseudogene/align/hg38_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/hg38_ago.txt \
  -bed $bedDir/formattedBed/hg38 \
  -meta $bedDir/hg38_all_dataset_info_table.txt \
  -output $baseDir/pseudogene/align/withSeq/hg38_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRpseudogenewithCLIP.hg38.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/sncRNA/align/hg38_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/hg38_ago.txt \
  -bed $bedDir/formattedBed/hg38 \
  -meta $bedDir/hg38_all_dataset_info_table.txt \
  -output $baseDir/sncRNA/align/withSeq/hg38_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRsncRNAwithCLIP.hg38.log 2>&1

