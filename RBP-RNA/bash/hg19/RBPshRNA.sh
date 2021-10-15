#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
logDir=$scripts/bash/RBP/log

#hg19
output=$base/RBPpredict/hg19

cd $output


# connected with shRBP data
RBPRef=$base/clip/rbsSeekerBed/reference/hg19_RBP_ENSEMBL.txt
shRBPFolder=$base/shRNA_RBPs_tsv/hg38

$scripts/RBPclusterAnnoShRBP.py -input anno/mRNA/hg19_RBPmRNA_gene.txt \
  -ref $RBPRef -shFolder $shRBPFolder -output anno/mRNA/hg19_RBPmRNA_gene.shRBP.txt

$scripts/RBPclusterAnnoShRBP.py -input anno/lncRNA/hg19_RBPlncRNA_gene.txt \
  -ref $RBPRef -shFolder $shRBPFolder -output anno/lncRNA/hg19_RBPlncRNA_gene.shRBP.txt

$scripts/RBPclusterAnnoShRBP.py -input anno/pseudogene/hg19_RBPpseudogene_gene.txt \
  -ref $RBPRef -shFolder $shRBPFolder -output anno/pseudogene/hg19_RBPpseudogene_gene.shRBP.txt

$scripts/RBPclusterAnnoShRBP.py -input anno/sncRNA/hg19_RBPsncRNA_gene.txt \
  -ref $RBPRef -shFolder $shRBPFolder -output anno/sncRNA/hg19_RBPsncRNA_gene.shRBP.txt

