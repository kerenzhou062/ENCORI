#!/bin/sh
bashDir=/data/liushun/starbase3/scripts/bash
logDir=$bashDir/miRother/log
cd ./mm10
starttime=`date +'%Y-%m-%d %H:%M:%S'`
#echo "filter miRotherPredict file"

./miRotherPredictFilter.sh > $logDir/miRotherPredictFilter.mm10.log 2>&1

echo "sorted miRotherPredict file intersect with AGO-CLIP data"
./miRotherWithCLIP.sh

echo "sorted miRotherPredict file intersect with degradome-seq data"

./miRotherWithDegradome.sh

echo "Integrate miRotherPredicts , AGO-CLIP and degradome-seq data"

./miRotherWithSeq.sh

echo "Generating microRNA-gene annotation"
./miRotherAnno.sh

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
runtime=$((end_seconds-start_seconds));
echo "Total runtime：${runtime}s"

