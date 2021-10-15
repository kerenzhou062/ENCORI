#!/bin/sh
bashDir=/data/liushun/starbase3/scripts/bash
logDir=$bashDir/miRmRNA/log
cd ./hg38
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "sort miRmRNApredict file"

./miRmRNApredictSort.sh > $logDir/miRmRNApredictSort.hg38.log 2>&1

echo "sorted miRmRNApredict file intersect with AGO-CLIP data"
./miRmRNAwithCLIP.sh

echo "sorted miRmRNApredict file intersect with degradome-seq data"

./miRmRNAwithDegradome.sh

echo "Integrate miRmRNApredicts , AGO-CLIP and degradome-seq data"

./miRmRNAwithSeq.sh

echo "Generating microRNA-cluster"
./miRmRNAcluster.sh

echo "Generating microRNA-pancancer"
./miRmRNApancancer.sh

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
runtime=$((end_seconds-start_seconds));
echo "Total runtime：${runtime}s"

