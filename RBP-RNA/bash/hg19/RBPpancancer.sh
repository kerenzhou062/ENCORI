#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
logDir=$scripts/bash/RBP/log
output=$base/RBPpredict/hg19/anno

RBPRef=$base/clip/rbsSeekerBed/reference/hg19_RBP_ENSEMBL.txt
RNAseq=$base/TCGA/common/RNAseq

#mRNA
cd $output/mRNA
awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2,$3)}}' hg19_RBPmRNA_gene.shRBP.txt | sort | uniq > hg19_RBP_geneID.txt
$scripts/RBPpanCancer.py -input hg19_RBP_geneID.txt -ref $RBPRef \
  -cpu 30 -cutoff 0.01 -seq $RNAseq -output hg19_RBP_pancancer.txt > $logDir/RBPpanCancer.mRNA.hg19.log 2>&1

$scripts/RBPannoPancancer.py -input hg19_RBPmRNA_gene.shRBP.txt \
  -pan hg19_RBP_pancancer.txt -output hg19_RBPmRNA.txt > $logDir/RBPannoPancancer.mRNA.hg19.log 2>&1

#lncRNA
cd $output/lncRNA
awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2,$3)}}' hg19_RBPlncRNA_gene.shRBP.txt | sort | uniq > hg19_RBP_geneID.txt
$scripts/RBPpanCancer.py -input hg19_RBP_geneID.txt -ref $RBPRef \
  -cpu 30 -cutoff 0.01 -seq $RNAseq -output hg19_RBP_pancancer.txt > $logDir/RBPpanCancer.lncRNA.hg19.log 2>&1

$scripts/RBPannoPancancer.py -input hg19_RBPlncRNA_gene.shRBP.txt \
  -pan hg19_RBP_pancancer.txt -output hg19_RBPlncRNA.txt > $logDir/RBPannoPancancer.lncRNA.hg19.log 2>&1

#pseudogene
cd $output/pseudogene
awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2,$3)}}' hg19_RBPpseudogene_gene.shRBP.txt | sort | uniq > hg19_RBP_geneID.txt
$scripts/RBPpanCancer.py -input hg19_RBP_geneID.txt -ref $RBPRef \
  -cpu 30 -cutoff 0.01 -seq $RNAseq -output hg19_RBP_pancancer.txt > $logDir/RBPpanCancer.pseudogene.hg19.log 2>&1

$scripts/RBPannoPancancer.py -input hg19_RBPpseudogene_gene.shRBP.txt \
  -pan hg19_RBP_pancancer.txt -output hg19_RBPpseudogene.txt > $logDir/RBPannoPancancer.pseudogene.hg19.log 2>&1

#sncRNA
cd $output/sncRNA
awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2,$3)}}' hg19_RBPsncRNA_gene.shRBP.txt | sort | uniq > hg19_RBP_geneID.txt
$scripts/RBPpanCancer.py -input hg19_RBP_geneID.txt -ref $RBPRef \
  -cpu 30 -cutoff 0.01 -seq $RNAseq -output hg19_RBP_pancancer.txt > $logDir/RBPpanCancer.sncRNA.hg19.log 2>&1

$scripts/RBPannoPancancer.py -input hg19_RBPsncRNA_gene.shRBP.txt \
  -pan hg19_RBP_pancancer.txt -output hg19_RBPsncRNA.txt > $logDir/RBPannoPancancer.sncRNA.hg19.log 2>&1
