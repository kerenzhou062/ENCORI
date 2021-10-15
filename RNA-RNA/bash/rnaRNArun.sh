#!/bin/sh
base=/data/liushun/starbase3
rriScanResDir="$base/rriScan/rriScan_new/filter"
reference="$base/reference"
scripts=$base/scripts
log=$base/scripts/bash/log

cd $base/rriScan

## hg38
#$scripts/rnaRNA.py -input $rriScanResDir/hg38_rnaRNApair.filter.txt \
#  -anno $reference/hg38/hg38.transcript.bed12 -ref hg38_rnaRNA_ref.txt -prefix hg38 -output ./

anno=$base/reference/hg38/typeClassAnno
$scripts/rnaRNAtypeAnno.py -input hg38_rnaRNA.txt \
  -anno $anno/hg38.transcript.mRNA.bed12 \
   $anno/hg38.transcript.lncRNA.bed12 \
   $anno/hg38.transcript.pseudogene.bed12 \
   $anno/hg38.transcript.sncRNA.bed12 \
   -output hg38_rnaRNA_gene.txt

## mm10
#$scripts/rnaRNA.py -input $rriScanResDir/mm10_rnaRNApair.filter.txt \
#  -anno $reference/mm10/mm10.transcript.bed12 -ref mm10_rnaRNA_ref.txt  -prefix mm10 -output ./

anno=$base/reference/mm10/typeClassAnno
$scripts/rnaRNAtypeAnno.py -input mm10_rnaRNA.txt \
  -anno $anno/mm10.transcript.mRNA.bed12 \
   $anno/mm10.transcript.lncRNA.bed12 \
   $anno/mm10.transcript.pseudogene.bed12 \
   $anno/mm10.transcript.sncRNA.bed12 \
   -output mm10_rnaRNA_gene.txt

echo "done"
