#!/bin/sh
bashDir=/data/liushun/starbase3/scripts/bash
logDir=$bashDir/RBP/log
cd $bashDir/RBP/hg38
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "generating RBP-cluster and annotation"
./RBPcluster.sh

echo "intersecting RBP-cluster-anno and shRBP data"
./RBPshRNA.sh

echo "link RBP-cluster-anno with pancancer"

./RBPpancancer.sh

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
runtime=$((end_seconds-start_seconds));
echo "Total runtime：${runtime}s"

