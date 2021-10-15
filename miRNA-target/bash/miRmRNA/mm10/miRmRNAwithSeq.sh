#!/bin/sh
base=/data/liushun/starbase3
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/mRNA/align/withSeq
logDir=$scripts/bash/miRmRNA/log

$scripts/miRmRNAwithSeq.py -clip $baseDir/mm10_all_sorted_cons_sb3_align.CLIP.txt \
  -degradome $baseDir/mm10_all_sorted_cons_sb3_align.degradome.txt \
  -output $baseDir/mm10_all_sorted_cons_sb3_align.withSeq.txt > $logDir/miRmRNAwithSeq.mm10.log 2>&1


awk 'BEGIN{OFS="\t";FS="\t";}
  {
    if($13>0){
      print $0;
    }
  }' $baseDir/mm10_all_sorted_cons_sb3_align.withSeq.txt > $baseDir/mm10_all_sorted_cons_sb3_align.withSeq.withCLIP.txt

awk 'BEGIN{OFS="\t";FS="\t";count=1;}
  {
    if(FNR==1){
      print "lineID","chromosome","start","end","bindID","merClass","strand","miRNAseq","align","targetSeq","tdmdScore","clipExpNum",
        "clipSiteNum","clipSiteID","degraExpNum","degraSiteNum","degraSiteID";
      print count,$1,$2,$3,$4,$5,$6,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18;
      count ++;
    }else{
      print count,$1,$2,$3,$4,$5,$6,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18;
      count ++;
    }
  }' $baseDir/mm10_all_sorted_cons_sb3_align.withSeq.withCLIP.txt > $baseDir/../../mm10_miRmRNA_align.txt

cd /data/liushun/starbase3/miRpredict/mRNA

$scripts/miRtargetCons.py -input mm10_miRmRNA_align.txt \
  -conscore $base/annotations/mm10/mm10.60way.phastCons.bedgraph \
  -output mm10_miRmRNA_align.cons.txt

mv mm10_miRmRNA_align.cons.txt mm10_miRmRNA_align.txt
