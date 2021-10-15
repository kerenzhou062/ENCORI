#!/bin/sh
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict
meta=/data/liushun/starbase3/clip/rbsSeekerBed/hg19_all_dataset_info_table.txt
logDir=$scripts/bash/miRother/log

for i in $baseDir/*;do
    if [[ -d $i ]];then
        baseName=${i##*/}
        if [[ $baseName != 'mRNA' ]]; then
            $scripts/miRotherWithSeq.py -clip $i/align/withSeq/hg19_cons_strict_site_sb3_id.CLIP.txt \
            -degradome $i/align/withSeq/hg19_cons_strict_site_sb3_id.degradome.txt \
            -meta $meta \
            -output $i/hg19_cons_strict_site_sb3_id.txt > $logDir/miR${baseName}withSeq.hg19.log 2>&1
        fi
    fi
done

