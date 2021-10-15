#!/bin/sh

cd /data/liushun/starbase3/scripts/bash/miRother
nohup ./pipeline-miRother.hg19.sh > ./log/pipeline-miRother.hg19.log 2>&1 &
nohup ./pipeline-miRother.hg38.sh > ./log/pipeline-miRother.hg38.log 2>&1 &
nohup ./pipeline-miRother.mm10.sh > ./log/pipeline-miRother.mm10.log 2>&1 &
