#!/bin/sh
base=/data/liushun/starbase3
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/mRNA/align/withSeq
logDir=$scripts/bash/miRmRNA/log

$scripts/miRmRNAwithSeq.py -clip $baseDir/hg19_all_sorted_cons_sb3_align.CLIP.txt \
  -degradome $baseDir/hg19_all_sorted_cons_sb3_align.degradome.txt \
  -output $baseDir/hg19_all_sorted_cons_sb3_align.withSeq.txt > $logDir/miRmRNAwithSeq.hg19.log 2>&1

awk 'BEGIN{OFS="\t";FS="\t";}
  {
    if($13>0){
      print $0;
    }
  }' $baseDir/hg19_all_sorted_cons_sb3_align.withSeq.txt > $baseDir/hg19_all_sorted_cons_sb3_align.withSeq.withCLIP.txt

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
  }' $baseDir/hg19_all_sorted_cons_sb3_align.withSeq.withCLIP.txt > $baseDir/../../hg19_miRmRNA_align.txt

cd /data/liushun/starbase3/miRpredict/mRNA

$scripts/miRtargetCons.py -input hg19_miRmRNA_align.txt \
  -conscore $base/annotations/hg19/hg19.100way.phastCons.bedgraph \
  -output hg19_miRmRNA_align.cons.txt

mv hg19_miRmRNA_align.cons.txt hg19_miRmRNA_align.txt
