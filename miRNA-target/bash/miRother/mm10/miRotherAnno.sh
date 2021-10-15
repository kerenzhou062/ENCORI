#!/bin/sh
base=/data/liushun/starbase3
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict
logDir=$scripts/bash/miRother/log
annoDir=/data/liushun/starbase3/reference/mm10/typeClassAnno

typeARR=( 'circRNA' 'lncRNA' 'sncRNA' 'pseudogene' )
for i in "${typeARR[@]}"*;do
    i="$baseDir/$i"
    if [[ -d $i ]];then
        awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR > 1){if($13 > 0){print $0}}}' $i/mm10_cons_strict_site_sb3_id.txt \
            > $i/mm10_cons_strict_site_sb3_id.withCLIP.txt;
    fi
done

cd $baseDir/circRNA
anno=$annoDir/mm10.transcript.circRNA.bed12
bedtools intersect -a mm10_cons_strict_site_sb3_id.withCLIP.txt -b $anno -s -split -wa -wb > mm10_cons_strict_site_sb3_id.withCLIP.anno.temp.txt
$scripts/miRotherAnno.py -circRNA -type bed12 -input mm10_cons_strict_site_sb3_id.withCLIP.anno.temp.txt -output mm10_miRcircRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' mm10_miRcircRNA.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > mm10_miRcircRNA_miRNA.txt

$scripts/miRtargetCons.py -input mm10_miRcircRNA.txt \
  -conscore $base/annotations/mm10/mm10.60way.phastCons.bedgraph \
  -output mm10_miRcircRNA.cons.txt

mv mm10_miRcircRNA.cons.txt mm10_miRcircRNA.txt


cd $baseDir/lncRNA
anno=$annoDir/mm10.transcript.lncRNA.bed12
bedtools intersect -a mm10_cons_strict_site_sb3_id.withCLIP.txt -b $anno -s -split -wa -wb > mm10_cons_strict_site_sb3_id.withCLIP.anno.temp.txt
$scripts/miRotherAnno.py -type bed12 -input mm10_cons_strict_site_sb3_id.withCLIP.anno.temp.txt -output mm10_miRlncRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' mm10_miRlncRNA.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > mm10_miRlncRNA_miRNA.txt

$scripts/miRtargetCons.py -input mm10_miRlncRNA.txt \
  -conscore $base/annotations/mm10/mm10.60way.phastCons.bedgraph \
  -output mm10_miRlncRNA.cons.txt

mv mm10_miRlncRNA.cons.txt mm10_miRlncRNA.txt


cd $baseDir/pseudogene
anno=$annoDir/mm10.transcript.pseudogene.bed12
bedtools intersect -a mm10_cons_strict_site_sb3_id.withCLIP.txt -b $anno -s -split -wa -wb > mm10_cons_strict_site_sb3_id.withCLIP.anno.temp.txt
$scripts/miRotherAnno.py -type bed12 -input mm10_cons_strict_site_sb3_id.withCLIP.anno.temp.txt -output mm10_miRpseudogene.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' mm10_miRpseudogene.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > mm10_miRpseudogene_miRNA.txt

$scripts/miRtargetCons.py -input mm10_miRpseudogene.txt \
  -conscore $base/annotations/mm10/mm10.60way.phastCons.bedgraph \
  -output mm10_miRpseudogene.cons.txt

mv mm10_miRpseudogene.cons.txt mm10_miRpseudogene.txt


cd $baseDir/sncRNA
anno=$annoDir/mm10.gene.sncRNA.bed
bedtools intersect -a mm10_cons_strict_site_sb3_id.withCLIP.txt -b $anno -s -split -wa -wb > mm10_cons_strict_site_sb3_id.withCLIP.anno.temp.txt
$scripts/miRotherAnno.py -type bed6 -input mm10_cons_strict_site_sb3_id.withCLIP.anno.temp.txt -output mm10_miRsncRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' mm10_miRsncRNA.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > mm10_miRsncRNA_miRNA.txt

$scripts/miRtargetCons.py -input mm10_miRsncRNA.txt \
  -conscore $base/annotations/mm10/mm10.60way.phastCons.bedgraph \
  -output mm10_miRsncRNA.cons.txt

mv mm10_miRsncRNA.cons.txt mm10_miRsncRNA.txt
