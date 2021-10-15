#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
logDir=$scripts/bash/RBP/log

#hg38
output=$base/RBPpredict/hg38

cd $output


# connected with shRBP data
RBPRef=$base/clip/rbsSeekerBed/reference/hg38_RBP_ENSEMBL.txt
shRBPFolder=$base/shRNA_RBPs_tsv/hg38

$scripts/RBPclusterAnnoShRBP.py -input anno/mRNA/hg38_RBPmRNA_gene.txt \
  -ref $RBPRef -shFolder $shRBPFolder -output anno/mRNA/hg38_RBPmRNA_gene.shRBP.txt

$scripts/RBPclusterAnnoShRBP.py -input anno/lncRNA/hg38_RBPlncRNA_gene.txt \
  -ref $RBPRef -shFolder $shRBPFolder -output anno/lncRNA/hg38_RBPlncRNA_gene.shRBP.txt

$scripts/RBPclusterAnnoShRBP.py -input anno/pseudogene/hg38_RBPpseudogene_gene.txt \
  -ref $RBPRef -shFolder $shRBPFolder -output anno/pseudogene/hg38_RBPpseudogene_gene.shRBP.txt

$scripts/RBPclusterAnnoShRBP.py -input anno/sncRNA/hg38_RBPsncRNA_gene.txt \
  -ref $RBPRef -shFolder $shRBPFolder -output anno/sncRNA/hg38_RBPsncRNA_gene.shRBP.txt

