#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
logDir=$scripts/bash/RBP/log

#ce10
bedDir=$base/clip/rbsSeekerBed/formattedBed/ce10
meta=$base/clip/rbsSeekerBed/ce10_all_dataset_info_table.txt
agoRef=$base/clip/rbsSeekerBed/reference/ce10_ago.txt
output=$base/RBPpredict/ce10

cd $output

rm -f catBed/*

# cat RBP
$scripts/RBPcatBed.py -cpu 10  -input $bedDir -agoRef $agoRef \
  -meta $meta -output catBed > $logDir/RBPcatSort.ce10.log 2>&1

# cluster RBP
rm -f clusterBed/*
$scripts/RBPclusterBed.py -cpu 10 -distance 20 -input catBed \
  -prefix CH -output clusterBed > $logDir/RBPclusterBed.ce10.log 2>&1

cat clusterBed/*.txt > ce10_RBP_cluster.txt

find ./anno -type f -name "*.txt" | xargs -I {} rm -f {}

# annotate RBP-mRNA
cd $output/anno/mRNA
rm -f ./*
anno=$base/reference/ce10/typeClassAnno/ce10.transcript.mRNA.bed12

bedtools intersect -a $output/ce10_RBP_cluster.txt -b $anno -s -split -wa -wb > ce10_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input ce10_RBP_cluster.anno.temp.txt \
  -type bed12 -cds -output ce10_RBPmRNA_cluster.txt > $logDir/RBPclusterAnno.ce10.mRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input ce10_RBPmRNA_cluster.txt -output ce10_RBPmRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' ce10_RBPmRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > ce10_RBPmRNA_RBP.txt

# annotate RBP-lncRNA
cd $output/anno/lncRNA
rm -f ./*
anno=$base/reference/ce10/typeClassAnno/ce10.transcript.lncRNA.bed12

bedtools intersect -a $output/ce10_RBP_cluster.txt -b $anno -s -split -wa -wb > ce10_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input ce10_RBP_cluster.anno.temp.txt \
  -type bed12 -output ce10_RBPlncRNA_cluster.txt > $logDir/RBPclusterAnno.ce10.lncRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input ce10_RBPlncRNA_cluster.txt -output ce10_RBPlncRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' ce10_RBPlncRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > ce10_RBPlncRNA_RBP.txt

# annotate RBP-pseudogene
cd $output/anno/pseudogene
rm -f ./*
anno=$base/reference/ce10/typeClassAnno/ce10.transcript.pseudogene.bed12

bedtools intersect -a $output/ce10_RBP_cluster.txt -b $anno -s -split -wa -wb > ce10_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input ce10_RBP_cluster.anno.temp.txt \
  -type bed12 -output ce10_RBPpseudogene_cluster.txt > $logDir/RBPclusterAnno.ce10.pseudogene.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input ce10_RBPpseudogene_cluster.txt -output ce10_RBPpseudogene.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' ce10_RBPpseudogene.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > ce10_RBPpseudogene_RBP.txt

# annotate RBP-sncRNA
cd $output/anno/sncRNA
rm -f ./*
anno=$base/reference/ce10/typeClassAnno/ce10.gene.sncRNA.bed

bedtools intersect -a $output/ce10_RBP_cluster.txt -b $anno -s -split -wa -wb > ce10_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input ce10_RBP_cluster.anno.temp.txt \
  -type bed6 -output ce10_RBPsncRNA_cluster.txt > $logDir/RBPclusterAnno.ce10.sncRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input ce10_RBPsncRNA_cluster.txt -output ce10_RBPsncRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' ce10_RBPsncRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > ce10_RBPsncRNA_RBP.txt

# annotate RBP-circRNA
cd $output/anno/circRNA
rm -f ./*
anno=$base/reference/ce10/typeClassAnno/ce10.transcript.circRNA.bed12

bedtools intersect -a $output/ce10_RBP_cluster.txt -b $anno -s -split -wa -wb > ce10_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input ce10_RBP_cluster.anno.temp.txt \
  -type bed12 -cds -circRNA -output ce10_RBPcircRNA_cluster.txt > $logDir/RBPclusterAnno.ce10.circRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input ce10_RBPcircRNA_cluster.txt -output ce10_RBPcircRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' ce10_RBPcircRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > ce10_RBPcircRNA_RBP.txt
