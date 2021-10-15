#!/bin/sh
base="/data/liushun/starbase3"
scripts=$base/scripts
baseDir=$base/miRpredict
bedDir=$base/clip/rbsSeekerBed
logDir=$scripts/bash/miRother/log

#circRNA,lncRNA,pseudogene,sncRNA
$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/lncRNA/align/mm10_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/mm10_ago.txt \
  -bed $bedDir/formattedBed/mm10 \
  -meta $bedDir/mm10_all_dataset_info_table.txt \
  -output $baseDir/lncRNA/align/withSeq/mm10_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRlncRNAwithCLIP.mm10.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/circRNA/align/mm10_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/mm10_ago.txt \
  -bed $bedDir/formattedBed/mm10 \
  -meta $bedDir/mm10_all_dataset_info_table.txt \
  -output $baseDir/circRNA/align/withSeq/mm10_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRcircRNAwithCLIP.mm10.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/pseudogene/align/mm10_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/mm10_ago.txt \
  -bed $bedDir/formattedBed/mm10 \
  -meta $bedDir/mm10_all_dataset_info_table.txt \
  -output $baseDir/pseudogene/align/withSeq/mm10_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRpseudogenewithCLIP.mm10.log 2>&1

$scripts/miRmRNAwithCLIP.py -cpu 30 -input $baseDir/sncRNA/align/mm10_cons_strict_site_sb3_id.bed \
  -agoRef $bedDir/reference/mm10_ago.txt \
  -bed $bedDir/formattedBed/mm10 \
  -meta $bedDir/mm10_all_dataset_info_table.txt \
  -output $baseDir/sncRNA/align/withSeq/mm10_cons_strict_site_sb3_id.CLIP.txt > $logDir/miRsncRNAwithCLIP.mm10.log 2>&1

