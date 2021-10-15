#!/bin/sh
base=/data/liushun/starbase3
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict
logDir=$scripts/bash/miRother/log
annoDir=/data/liushun/starbase3/reference/hg38/typeClassAnno

typeARR=( 'circRNA' 'lncRNA' 'sncRNA' 'pseudogene' )
for i in "${typeARR[@]}"*;do
    i="$baseDir/$i"
    if [[ -d $i ]];then
        awk 'BEGIN{FS="\t";OFS="\t";}{if(FNR > 1){if($13 > 0){print $0}}}' $i/hg38_cons_strict_site_sb3_id.txt \
            > $i/hg38_cons_strict_site_sb3_id.withCLIP.txt;
    fi
done

cd $baseDir/lncRNA
anno=$annoDir/hg38.transcript.lncRNA.bed12
bedtools intersect -a hg38_cons_strict_site_sb3_id.withCLIP.txt -b $anno -s -split -wa -wb > hg38_cons_strict_site_sb3_id.withCLIP.anno.temp.txt
$scripts/miRotherAnno.py -type bed12 -input hg38_cons_strict_site_sb3_id.withCLIP.anno.temp.txt -output hg38_cons_strict_site_sb3_id.withCLIP.anno.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' hg38_cons_strict_site_sb3_id.withCLIP.anno.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > hg38_miRlncRNA_miRNA.txt

cd $baseDir/pseudogene
anno=$annoDir/hg38.transcript.pseudogene.bed12
bedtools intersect -a hg38_cons_strict_site_sb3_id.withCLIP.txt -b $anno -s -split -wa -wb > hg38_cons_strict_site_sb3_id.withCLIP.anno.temp.txt
$scripts/miRotherAnno.py -type bed12 -input hg38_cons_strict_site_sb3_id.withCLIP.anno.temp.txt -output hg38_cons_strict_site_sb3_id.withCLIP.anno.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' hg38_cons_strict_site_sb3_id.withCLIP.anno.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > hg38_miRpseudogene_miRNA.txt

cd $baseDir/sncRNA
anno=$annoDir/hg38.gene.sncRNA.bed
bedtools intersect -a hg38_cons_strict_site_sb3_id.withCLIP.txt -b $anno -s -split -wa -wb > hg38_cons_strict_site_sb3_id.withCLIP.anno.temp.txt
$scripts/miRotherAnno.py -type bed6 -input hg38_cons_strict_site_sb3_id.withCLIP.anno.temp.txt -output hg38_cons_strict_site_sb3_id.withCLIP.anno.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' hg38_cons_strict_site_sb3_id.withCLIP.anno.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > hg38_miRsncRNA_miRNA.txt

cd $baseDir/circRNA
anno=$annoDir/hg38.transcript.circRNA.bed12
bedtools intersect -a hg38_cons_strict_site_sb3_id.withCLIP.txt -b $anno -s -split -wa -wb > hg38_cons_strict_site_sb3_id.withCLIP.anno.temp.txt
$scripts/miRotherAnno.py -circRNA -type bed12 -input hg38_cons_strict_site_sb3_id.withCLIP.anno.temp.txt -output hg38_miRcircRNA.txt

awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR>1){print($8,$9)}}' hg38_miRcircRNA.txt | sort | uniq |\
  awk 'BEGIN{OFS="\t";FS="\t";}{if(FNR==1){print "miRNAid", "miRNAname";print $0;}else{print $0}}' > hg38_miRcircRNA_miRNA.txt

$scripts/miRtargetCons.py -input hg38_miRcircRNA.txt \
  -conscore $base/annotations/hg38/hg38.phastCons100way.bedgraph \
  -output hg38_miRcircRNA.cons.txt

mv hg38_miRcircRNA.cons.txt hg38_miRcircRNA.txt
