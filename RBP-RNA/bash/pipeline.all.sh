#!/bin/sh

cd /data/liushun/starbase3/scripts/bash/RBP
nohup ./pipeline-RBP.hg19.sh > ./log/pipeline-RBP.hg19.log 2>&1 &
nohup ./pipeline-RBP.hg38.sh > ./log/pipeline-RBP.hg38.log 2>&1 &
nohup ./pipeline-RBP.mm10.sh > ./log/pipeline-RBP.mm10.log 2>&1 &
nohup ./pipeline-RBP.ce10.sh > ./log/pipeline-RBP.ce10.log 2>&1 &
nohup ./pipeline-RBP.dm6.sh > ./log/pipeline-RBP.dm6.log 2>&1 &
nohup ./pipeline-RBP.sacCer3.sh > ./log/pipeline-RBP.sacCer3.log 2>&1 &
