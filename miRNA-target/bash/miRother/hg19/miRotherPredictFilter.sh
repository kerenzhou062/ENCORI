#!/bin/sh
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict
fasta=/data/liushun/starbase3/reference/hg19/hg19.fa
primaryMirna=/data/liushun/starbase3/reference/hg19/miRBase.v20.primary.bed
starttime=`date +'%Y-%m-%d %H:%M:%S'`

# generate reference file
echo "filtering original predicting data"
typeARR=( 'circRNA' 'lncRNA' 'sncRNA' 'pseudogene' )
for i in "${typeARR[@]}"*;do
    i="$baseDir/$i"
    if [[ -d $i ]];then
        baseName=${i##*/}
        for j in $i/align/original/hg19_*.bed;do
            fileName=${j##*/}
            $scripts/miRtargetFilter.py -input $j -fasta $fasta -output $i/align/$fileName.tmp;
            bedtools intersect -a $i/align/$fileName.tmp -b $primaryMirna -s -v > $i/align/$fileName.tmp2
            bedtools intersect -a $primaryMirna -b $i/align/$fileName.tmp -s -wa -wb | \
              awk '{IGNORECASE=1;split($4, arr, ":"); if ($14 !~ arr[2]) {print}}' | cut -d $'\t' --complement -f-6 > $i/align/$fileName.tmp3
            cat $i/align/$fileName.tmp2 $i/align/$fileName.tmp3 > $i/align/$fileName
            rm -f $i/align/$fileName.tmp $i/align/$fileName.tmp2 $i/align/$fileName.tmp3
        done;
    fi
done

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
runtime=$((end_seconds-start_seconds));
echo "Total runtime：${runtime}s"

