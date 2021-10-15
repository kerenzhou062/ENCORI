#!/bin/sh

cd /data/liushun/starbase3/scripts/bash/miRmRNA
nohup ./pipeline-miRmRNA.hg19.sh > ./log/pipeline-miRmRNA.hg19.log 2>&1 &
nohup ./pipeline-miRmRNA.hg38.sh > ./log/pipeline-miRmRNA.hg38.log 2>&1 &
nohup ./pipeline-miRmRNA.mm10.sh > ./log/pipeline-miRmRNA.mm10.log 2>&1 &
