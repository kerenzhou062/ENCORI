#!/bin/sh
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/mRNA
logDir=$scripts/bash/miRmRNA/log
meta=/data/liushun/starbase3/clip/rbsSeekerBed/hg19_all_dataset_info_table.txt
anno=/data/liushun/starbase3/reference/hg19/typeClassAnno/hg19.transcript.mRNA.bed12

cd $baseDir

$scripts/miRmRNAcluster.py -input $baseDir/align/withSeq/hg19_all_sorted_cons_sb3_align.withSeq.txt \
  -ref $baseDir/reference/hg19_bindingID_reference.txt \
  -prefix CH \
  -meta $meta \
  -output hg19_clusterID_cons_sb3.txt > $logDir/miRmRNAcluster.hg19.log 2>&1

$scripts/miRmRNAcluster.py -input $baseDir/align/withSeq/hg19_all_sorted_cons_sb3_align.withSeq.withCLIP.txt \
  -ref $baseDir/reference/hg19_bindingID_reference.txt \
  -prefix CH \
  -meta $meta \
  -output hg19_clusterID_cons_sb3.withCLIP.txt > $logDir/miRmRNAcluster.hg19.log 2>&1

awk 'BEGEIN{OFS="\t";FS="\t";}{if(FNR>1){print $0}}' hg19_clusterID_cons_sb3.withCLIP.txt | \
  bedtools intersect -a stdin -b $anno -s -split -wa -wb > hg19_clusterID_cons_sb3.withCLIP.anno.temp.txt

$scripts/miRmRNAclusterAnno.py -input hg19_clusterID_cons_sb3.withCLIP.anno.temp.txt -output hg19_clusterID_cons_sb3.withCLIP.anno.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' hg19_clusterID_cons_sb3.withCLIP.anno.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > hg19_miRmRNA_miRNA.txt

