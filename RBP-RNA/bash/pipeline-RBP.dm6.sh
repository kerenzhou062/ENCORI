#!/bin/sh
bashDir=/data/liushun/starbase3/scripts/bash
logDir=$bashDir/RBP/log
cd $bashDir/RBP/dm6
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "generating RBP-cluster and annotation"
./RBPcluster.sh

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
runtime=$((end_seconds-start_seconds));
echo "Total runtime：${runtime}s"

