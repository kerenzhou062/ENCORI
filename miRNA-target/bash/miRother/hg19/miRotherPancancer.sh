#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
baseDir=$base/miRpredict
logDir=$scripts/bash/miRother/log

miRNAseq=$base/TCGA/common/miRNAseq
RNAseq=$base/TCGA/common/RNAseq

#lncRNA
cd $baseDir/lncRNA

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9,$24)}}' hg19_cons_strict_site_sb3_id.withCLIP.anno.txt | \
  sort | uniq > hg19_miRlncRNA_geneID.txt

$scripts/miRNApanCancer.py -input hg19_miRlncRNA_geneID.txt \
  -cpu 30 -cutoff 0.01 -miRNAseq $miRNAseq -RNAseq $RNAseq \
  -output hg19_lncRNApancancer.txt > $logDir/miRlncRNApancancer.hg19.log 2>&1

$scripts/miRannoPancancer.py -input hg19_cons_strict_site_sb3_id.withCLIP.anno.txt \
  -type lncRNA -pan hg19_lncRNApancancer.txt -output hg19_miRlncRNA.txt

$scripts/miRtargetCons.py -input hg19_miRlncRNA.txt \
  -conscore $base/annotations/hg19/hg19.100way.phastCons.bedgraph \
  -output hg19_miRlncRNA.cons.txt

mv hg19_miRlncRNA.cons.txt hg19_miRlncRNA.txt

#pseudogene
cd $baseDir/pseudogene

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9,$24)}}' hg19_cons_strict_site_sb3_id.withCLIP.anno.txt | \
  sort | uniq > hg19_miRpseudogene_geneID.txt

$scripts/miRNApanCancer.py -input hg19_miRpseudogene_geneID.txt \
  -cpu 30 -cutoff 0.01 -miRNAseq $miRNAseq -RNAseq $RNAseq \
  -output hg19_pseudogenepancancer.txt > $logDir/miRpseudogenepancancer.hg19.log 2>&1

$scripts/miRannoPancancer.py -input hg19_cons_strict_site_sb3_id.withCLIP.anno.txt \
  -type pseudogene -pan hg19_pseudogenepancancer.txt -output hg19_miRpseudogene.txt

$scripts/miRtargetCons.py -input hg19_miRpseudogene.txt \
  -conscore $base/annotations/hg19/hg19.100way.phastCons.bedgraph \
  -output hg19_miRpseudogene.cons.txt

mv hg19_miRpseudogene.cons.txt hg19_miRpseudogene.txt

#sncRNA
cd $baseDir/sncRNA

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9,$24)}}' hg19_cons_strict_site_sb3_id.withCLIP.anno.txt | \
  sort | uniq > hg19_miRsncRNA_geneID.txt

$scripts/miRNApanCancer.py -input hg19_miRsncRNA_geneID.txt \
  -cpu 30 -cutoff 0.01 -miRNAseq $miRNAseq -RNAseq $RNAseq \
  -output hg19_sncRNApancancer.txt > $logDir/miRsncRNApancancer.hg19.log 2>&1

$scripts/miRannoPancancer.py -input hg19_cons_strict_site_sb3_id.withCLIP.anno.txt \
  -type sncRNA -pan hg19_sncRNApancancer.txt -output hg19_miRsncRNA.txt

$scripts/miRtargetCons.py -input hg19_miRsncRNA.txt \
  -conscore $base/annotations/hg19/hg19.100way.phastCons.bedgraph \
  -output hg19_miRsncRNA.cons.txt

mv hg19_miRsncRNA.cons.txt hg19_miRsncRNA.txt
