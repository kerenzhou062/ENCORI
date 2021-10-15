#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
baseDir=$base/miRpredict/mRNA
logDir=$scripts/bash/miRmRNA/log

miRNAseq=$base/TCGA/common/miRNAseq
RNAseq=$base/TCGA/common/RNAseq

cd $baseDir

#awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9,$44)}}' hg19_clusterID_cons_sb3.withCLIP.anno.txt | sort | uniq > hg19_miRmRNA_geneID.txt
#
#$scripts/miRNApanCancer.py -input hg19_miRmRNA_geneID.txt \
#  -cpu 30 -cutoff 0.01 -miRNAseq $miRNAseq -RNAseq $RNAseq \
#  -output hg19_mRNApancancer.txt > $logDir/miRNApanCancer.hg19.log 2>&1

$scripts/miRannoPancancer.py -input hg19_clusterID_cons_sb3.withCLIP.anno.txt \
  -type mRNA -pan hg19_mRNApancancer.txt -output hg19_miRmRNA.txt

$scripts/miRtargetCons.py -input hg19_miRmRNA.txt \
  -conscore $base/annotations/hg19/hg19.100way.phastCons.bedgraph \
  -output hg19_miRmRNA.cons.txt

mv hg19_miRmRNA.cons.txt hg19_miRmRNA.txt
