#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
baseDir=$base/miRpredict
logDir=$scripts/bash/miRother/log

miRNAseq=$base/TCGA/common/miRNAseq
RNAseq=$base/TCGA/common/RNAseq

#lncRNA
cd $baseDir/lncRNA

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9,$24)}}' hg38_cons_strict_site_sb3_id.withCLIP.anno.txt | \
  sort | uniq > hg38_miRlncRNA_geneID.txt

$scripts/miRNApanCancer.py -input hg38_miRlncRNA_geneID.txt \
  -cpu 30 -cutoff 0.01 -miRNAseq $miRNAseq -RNAseq $RNAseq \
  -output hg38_lncRNApancancer.txt > $logDir/miRlncRNApancancer.hg38.log 2>&1

$scripts/miRannoPancancer.py -input hg38_cons_strict_site_sb3_id.withCLIP.anno.txt \
  -type lncRNA -pan hg38_lncRNApancancer.txt -output hg38_miRlncRNA.txt

$scripts/miRtargetCons.py -input hg38_miRlncRNA.txt \
  -conscore $base/annotations/hg38/hg38.phastCons100way.bedgraph \
  -output hg38_miRlncRNA.cons.txt

mv hg38_miRlncRNA.cons.txt hg38_miRlncRNA.txt

#pseudogene
cd $baseDir/pseudogene

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9,$24)}}' hg38_cons_strict_site_sb3_id.withCLIP.anno.txt | \
  sort | uniq > hg38_miRpseudogene_geneID.txt

$scripts/miRNApanCancer.py -input hg38_miRpseudogene_geneID.txt \
  -cpu 30 -cutoff 0.01 -miRNAseq $miRNAseq -RNAseq $RNAseq \
  -output hg38_pseudogenepancancer.txt > $logDir/miRpseudogenepancancer.hg38.log 2>&1

$scripts/miRannoPancancer.py -input hg38_cons_strict_site_sb3_id.withCLIP.anno.txt \
  -type pseudogene -pan hg38_pseudogenepancancer.txt -output hg38_miRpseudogene.txt

$scripts/miRtargetCons.py -input hg38_miRpseudogene.txt \
  -conscore $base/annotations/hg38/hg38.phastCons100way.bedgraph \
  -output hg38_miRpseudogene.cons.txt

mv hg38_miRpseudogene.cons.txt hg38_miRpseudogene.txt

#sncRNA
cd $baseDir/sncRNA

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9,$24)}}' hg38_cons_strict_site_sb3_id.withCLIP.anno.txt | \
  sort | uniq > hg38_miRsncRNA_geneID.txt

$scripts/miRNApanCancer.py -input hg38_miRsncRNA_geneID.txt \
  -cpu 30 -cutoff 0.01 -miRNAseq $miRNAseq -RNAseq $RNAseq \
  -output hg38_sncRNApancancer.txt > $logDir/miRsncRNApancancer.hg38.log 2>&1

$scripts/miRannoPancancer.py -input hg38_cons_strict_site_sb3_id.withCLIP.anno.txt \
  -type sncRNA -pan hg38_sncRNApancancer.txt -output hg38_miRsncRNA.txt

$scripts/miRtargetCons.py -input hg38_miRsncRNA.txt \
  -conscore $base/annotations/hg38/hg38.phastCons100way.bedgraph \
  -output hg38_miRsncRNA.cons.txt

mv hg38_miRsncRNA.cons.txt hg38_miRsncRNA.txt
