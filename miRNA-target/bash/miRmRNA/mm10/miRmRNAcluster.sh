#!/bin/sh
base=/data/liushun/starbase3/
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/mRNA
logDir=$scripts/bash/miRmRNA/log
meta=/data/liushun/starbase3/clip/rbsSeekerBed/mm10_all_dataset_info_table.txt
anno=/data/liushun/starbase3/reference/mm10/typeClassAnno/mm10.transcript.mRNA.bed12

cd $baseDir

$scripts/miRmRNAcluster.py -input $baseDir/align/withSeq/mm10_all_sorted_cons_sb3_align.withSeq.txt \
  -ref $baseDir/reference/mm10_bindingID_reference.txt \
  -prefix CM \
  -meta $meta \
  -output mm10_clusterID_cons_sb3.txt > $logDir/miRmRNAcluster.mm10.log 2>&1

$scripts/miRmRNAcluster.py -input $baseDir/align/withSeq/mm10_all_sorted_cons_sb3_align.withSeq.withCLIP.txt \
  -ref $baseDir/reference/mm10_bindingID_reference.txt \
  -prefix CM \
  -meta $meta \
  -output mm10_clusterID_cons_sb3.withCLIP.txt > $logDir/miRmRNAcluster.mm10.log 2>&1

awk 'BEGEIN{OFS="\t";FS="\t";}{if(FNR>1){print $0}}' mm10_clusterID_cons_sb3.withCLIP.txt | \
  bedtools intersect -a stdin -b $anno -s -split -wa -wb > mm10_clusterID_cons_sb3.withCLIP.anno.temp.txt

$scripts/miRmRNAclusterAnno.py -input mm10_clusterID_cons_sb3.withCLIP.anno.temp.txt -output mm10_miRmRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' mm10_miRmRNA.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > mm10_miRmRNA_miRNA.txt

$scripts/miRtargetCons.py -input mm10_miRmRNA.txt \
  -conscore $base/annotations/mm10/mm10.60way.phastCons.bedgraph \
  -output mm10_miRmRNA.cons.txt

mv mm10_miRmRNA.cons.txt mm10_miRmRNA.txt
