#!/bin/sh
base=/data/liushun/starbase3
scripts=$base/scripts
logDir=$scripts/bash/RBP/log

#hg19
bedDir=$base/clip/rbsSeekerBed/formattedBed/hg19
meta=$base/clip/rbsSeekerBed/hg19_all_dataset_info_table.txt
agoRef=$base/clip/rbsSeekerBed/reference/hg19_ago.txt
output=$base/RBPpredict/hg19

cd $output

rm -f catBed/*

# cat RBP
$scripts/RBPcatBed.py -cpu 10  -input $bedDir -agoRef $agoRef \
  -meta $meta -output catBed > $logDir/RBPcatSort.hg19.log 2>&1

# cluster RBP
rm -f clusterBed/*
$scripts/RBPclusterBed.py -cpu 10 -distance 20 -input catBed \
  -prefix CH -output clusterBed > $logDir/RBPclusterBed.hg19.log 2>&1

cat clusterBed/*.txt > hg19_RBP_cluster.txt

find ./anno -type f -name "*.txt" | xargs -I {} rm -f {}

# annotate RBP-mRNA
cd $output/anno/mRNA
rm -f ./*
anno=$base/reference/hg19/typeClassAnno/hg19.transcript.mRNA.bed12

bedtools intersect -a $output/hg19_RBP_cluster.txt -b $anno -s -split -wa -wb > hg19_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input hg19_RBP_cluster.anno.temp.txt \
  -type bed12 -cds -output hg19_RBPmRNA_cluster.txt > $logDir/RBPclusterAnno.hg19.mRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input hg19_RBPmRNA_cluster.txt -output hg19_RBPmRNA_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' hg19_RBPmRNA_gene.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > hg19_RBPmRNA_RBP.txt

# annotate RBP-lncRNA
cd $output/anno/lncRNA
rm -f ./*
anno=$base/reference/hg19/typeClassAnno/hg19.transcript.lncRNA.bed12

bedtools intersect -a $output/hg19_RBP_cluster.txt -b $anno -s -split -wa -wb > hg19_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input hg19_RBP_cluster.anno.temp.txt \
  -type bed12 -output hg19_RBPlncRNA_cluster.txt > $logDir/RBPclusterAnno.hg19.lncRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input hg19_RBPlncRNA_cluster.txt -output hg19_RBPlncRNA_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' hg19_RBPlncRNA_gene.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > hg19_RBPlncRNA_RBP.txt

# annotate RBP-pseudogene
cd $output/anno/pseudogene
rm -f ./*
anno=$base/reference/hg19/typeClassAnno/hg19.transcript.pseudogene.bed12

bedtools intersect -a $output/hg19_RBP_cluster.txt -b $anno -s -split -wa -wb > hg19_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input hg19_RBP_cluster.anno.temp.txt \
  -type bed12 -output hg19_RBPpseudogene_cluster.txt > $logDir/RBPclusterAnno.hg19.pseudogene.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input hg19_RBPpseudogene_cluster.txt -output hg19_RBPpseudogene_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' hg19_RBPpseudogene_gene.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > hg19_RBPpseudogene_RBP.txt

# annotate RBP-sncRNA
cd $output/anno/sncRNA
rm -f ./*
anno=$base/reference/hg19/typeClassAnno/hg19.gene.sncRNA.bed

bedtools intersect -a $output/hg19_RBP_cluster.txt -b $anno -s -split -wa -wb > hg19_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input hg19_RBP_cluster.anno.temp.txt \
  -type bed6 -output hg19_RBPsncRNA_cluster.txt > $logDir/RBPclusterAnno.hg19.sncRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input hg19_RBPsncRNA_cluster.txt -output hg19_RBPsncRNA_gene.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' hg19_RBPsncRNA_gene.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > hg19_RBPsncRNA_RBP.txt

# annotate RBP-circRNA
cd $output/anno/circRNA
rm -f ./*
anno=$base/reference/hg19/typeClassAnno/hg19.transcript.circRNA.bed12

bedtools intersect -a $output/hg19_RBP_cluster.txt -b $anno -s -split -wa -wb > hg19_RBP_cluster.anno.temp.txt

$scripts/RBPclusterBedAnno.py -input hg19_RBP_cluster.anno.temp.txt \
  -type bed12 -cds -circRNA -output hg19_RBPcircRNA_cluster.txt > $logDir/RBPclusterAnno.hg19.circRNA.log 2>&1

$scripts/RBPclusterBedGeneMerge.py -input hg19_RBPcircRNA_cluster.txt -output hg19_RBPcircRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR > 1){print($2)}}' hg19_RBPcircRNA.txt | sort | uniq | \
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "RBPname";print $0;}else{print $0}}' > hg19_RBPcircRNA_RBP.txt
