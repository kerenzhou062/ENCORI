#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
logDir=$scripts/bash/RBP/log

#sacCer3
bedDir=$base/clip/rbsSeekerBed/formattedBed/sacCer3
meta=$base/clip/rbsSeekerBed/sacCer3_all_dataset_info_table.txt
agoRef=$base/clip/rbsSeekerBed/reference/sacCer3_ago.txt
output=$base/RBPpredict/sacCer3

cd $output

rm -f catBed/*

# cat RBP
$scripts/RBPcatBed.py -cpu 10  -input $bedDir -agoRef $agoRef \
  -meta $meta -output catBed > $logDir/RBPcatSort.sacCer3.log 2>&1

# cluster RBP
rm -f clusterBed/*
$scripts/RBPclusterBed.py -cpu 10 -distance 20 -input catBed \
  -prefix CH -output clusterBed > $logDir/RBPclusterBed.sacCer3.log 2>&1

cat clusterBed/*.txt > sacCer3_RBP_cluster.txt

find ./anno -type f -name "*.txt" | xargs -I {} rm -f {}

# annotate RBP-mRNA
cd $output/anno/mRNA
rm -f ./*
anno=$base/reference/sacCer3/typeClassAnno/sacCer3.transcript.mRNA.bed12

bedtools intersect -a $output/sacCer3_RBP_cluster.txt -b $anno -s -split -wa -wb > sacCer3_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input sacCer3_RBP_cluster.anno.temp.txt \
  -type bed12 -cds -output sacCer3_RBPmRNA_cluster.txt > $logDir/RBPclusterAnno.sacCer3.mRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input sacCer3_RBPmRNA_cluster.txt -output sacCer3_RBPmRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' sacCer3_RBPmRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > sacCer3_RBPmRNA_RBP.txt

# annotate RBP-lncRNA
cd $output/anno/lncRNA
rm -f ./*
anno=$base/reference/sacCer3/typeClassAnno/sacCer3.transcript.lncRNA.bed12

bedtools intersect -a $output/sacCer3_RBP_cluster.txt -b $anno -s -split -wa -wb > sacCer3_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input sacCer3_RBP_cluster.anno.temp.txt \
  -type bed12 -output sacCer3_RBPlncRNA_cluster.txt > $logDir/RBPclusterAnno.sacCer3.lncRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input sacCer3_RBPlncRNA_cluster.txt -output sacCer3_RBPlncRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' sacCer3_RBPlncRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > sacCer3_RBPlncRNA_RBP.txt

# annotate RBP-pseudogene
cd $output/anno/pseudogene
rm -f ./*
anno=$base/reference/sacCer3/typeClassAnno/sacCer3.transcript.pseudogene.bed12

bedtools intersect -a $output/sacCer3_RBP_cluster.txt -b $anno -s -split -wa -wb > sacCer3_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input sacCer3_RBP_cluster.anno.temp.txt \
  -type bed12 -output sacCer3_RBPpseudogene_cluster.txt > $logDir/RBPclusterAnno.sacCer3.pseudogene.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input sacCer3_RBPpseudogene_cluster.txt -output sacCer3_RBPpseudogene.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' sacCer3_RBPpseudogene.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > sacCer3_RBPpseudogene_RBP.txt

# annotate RBP-sncRNA
cd $output/anno/sncRNA
rm -f ./*
anno=$base/reference/sacCer3/typeClassAnno/sacCer3.gene.sncRNA.bed

bedtools intersect -a $output/sacCer3_RBP_cluster.txt -b $anno -s -split -wa -wb > sacCer3_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input sacCer3_RBP_cluster.anno.temp.txt \
  -type bed6 -output sacCer3_RBPsncRNA_cluster.txt > $logDir/RBPclusterAnno.sacCer3.sncRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input sacCer3_RBPsncRNA_cluster.txt -output sacCer3_RBPsncRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' sacCer3_RBPsncRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > sacCer3_RBPsncRNA_RBP.txt
