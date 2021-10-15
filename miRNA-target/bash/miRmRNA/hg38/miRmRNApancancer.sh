#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
baseDir=$base/miRpredict/mRNA
logDir=$scripts/bash/miRmRNA/log

miRNAseq=$base/TCGA/common/miRNAseq
RNAseq=$base/TCGA/common/RNAseq

cd $baseDir

#awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9,$44)}}' hg38_clusterID_cons_sb3.withCLIP.anno.txt | sort | uniq > hg38_miRmRNA_geneID.txt
#
#$scripts/miRNApanCancer.py -input hg38_miRmRNA_geneID.txt \
#  -cpu 30 -cutoff 0.01 -miRNAseq $miRNAseq -RNAseq $RNAseq \
#  -output hg38_mRNApancancer.txt > $logDir/miRNApanCancer.hg38.log 2>&1

$scripts/miRannoPancancer.py -input hg38_clusterID_cons_sb3.withCLIP.anno.txt \
  -type mRNA -pan hg38_mRNApancancer.txt -output hg38_miRmRNA.txt

$scripts/miRtargetCons.py -input hg38_miRmRNA.txt \
  -conscore $base/annotations/hg38/hg38.phastCons100way.bedgraph \
  -output hg38_miRmRNA.cons.txt

mv hg38_miRmRNA.cons.txt hg38_miRmRNA.txt
