#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
logDir=$scripts/bash/RBP/log

#dm6
bedDir=$base/clip/rbsSeekerBed/formattedBed/dm6
meta=$base/clip/rbsSeekerBed/dm6_all_dataset_info_table.txt
agoRef=$base/clip/rbsSeekerBed/reference/dm6_ago.txt
output=$base/RBPpredict/dm6

cd $output

rm -f catBed/*

# cat RBP
$scripts/RBPcatBed.py -cpu 10  -input $bedDir -agoRef $agoRef \
  -meta $meta -output catBed > $logDir/RBPcatSort.dm6.log 2>&1

# cluster RBP
rm -f clusterBed/*
$scripts/RBPclusterBed.py -cpu 10 -distance 20 -input catBed \
  -prefix CH -output clusterBed > $logDir/RBPclusterBed.dm6.log 2>&1

cat clusterBed/*.txt > dm6_RBP_cluster.txt

find ./anno -type f -name "*.txt" | xargs -I {} rm -f {}

# annotate RBP-mRNA
cd $output/anno/mRNA
rm -f ./*
anno=$base/reference/dm6/typeClassAnno/dm6.transcript.mRNA.bed12

bedtools intersect -a $output/dm6_RBP_cluster.txt -b $anno -s -split -wa -wb > dm6_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input dm6_RBP_cluster.anno.temp.txt \
  -type bed12 -cds -output dm6_RBPmRNA_cluster.txt > $logDir/RBPclusterAnno.dm6.mRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input dm6_RBPmRNA_cluster.txt -output dm6_RBPmRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' dm6_RBPmRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > dm6_RBPmRNA_RBP.txt

# annotate RBP-lncRNA
cd $output/anno/lncRNA
rm -f ./*
anno=$base/reference/dm6/typeClassAnno/dm6.transcript.lncRNA.bed12

bedtools intersect -a $output/dm6_RBP_cluster.txt -b $anno -s -split -wa -wb > dm6_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input dm6_RBP_cluster.anno.temp.txt \
  -type bed12 -output dm6_RBPlncRNA_cluster.txt > $logDir/RBPclusterAnno.dm6.lncRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input dm6_RBPlncRNA_cluster.txt -output dm6_RBPlncRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' dm6_RBPlncRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > dm6_RBPlncRNA_RBP.txt

# annotate RBP-pseudogene
cd $output/anno/pseudogene
rm -f ./*
anno=$base/reference/dm6/typeClassAnno/dm6.transcript.pseudogene.bed12

bedtools intersect -a $output/dm6_RBP_cluster.txt -b $anno -s -split -wa -wb > dm6_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input dm6_RBP_cluster.anno.temp.txt \
  -type bed12 -output dm6_RBPpseudogene_cluster.txt > $logDir/RBPclusterAnno.dm6.pseudogene.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input dm6_RBPpseudogene_cluster.txt -output dm6_RBPpseudogene.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' dm6_RBPpseudogene.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > dm6_RBPpseudogene_RBP.txt

# annotate RBP-sncRNA
cd $output/anno/sncRNA
rm -f ./*
anno=$base/reference/dm6/typeClassAnno/dm6.gene.sncRNA.bed

bedtools intersect -a $output/dm6_RBP_cluster.txt -b $anno -s -split -wa -wb > dm6_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input dm6_RBP_cluster.anno.temp.txt \
  -type bed6 -output dm6_RBPsncRNA_cluster.txt > $logDir/RBPclusterAnno.dm6.sncRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input dm6_RBPsncRNA_cluster.txt -output dm6_RBPsncRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' dm6_RBPsncRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > dm6_RBPsncRNA_RBP.txt
