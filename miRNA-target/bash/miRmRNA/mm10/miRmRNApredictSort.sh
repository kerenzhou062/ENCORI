#!/bin/sh
scripts=/data/liushun/starbase3/scripts
baseDir=/data/liushun/starbase3/miRpredict/mRNA
fasta=/data/liushun/starbase3/reference/mm10/mm10.fa
primaryMirna=/data/liushun/starbase3/reference/mm10/miRBase.v22.primary.bed
starttime=`date +'%Y-%m-%d %H:%M:%S'`

# generate reference file
echo "filtering original predicting data"

for i in $baseDir/align/original/mm10_*.bed;do
    baseName=${i##*/}
    ## filter by type and primary miRNA gene locus
    $scripts/miRtargetFilter.py -input $i -fasta $fasta -output $baseDir/align/$baseName.tmp;
    bedtools intersect -a $baseDir/align/$baseName.tmp -b $primaryMirna -s -v > $baseDir/align/$baseName.tmp2
    bedtools intersect -a $primaryMirna -b $baseDir/align/$baseName.tmp -s -wa -wb | \
      awk '{IGNORECASE=1;split($4, arr, ":"); if ($14 !~ arr[2]) {print}}' | cut -d $'\t' --complement -f-6 > $baseDir/align/$baseName.tmp3
    cat $baseDir/align/$baseName.tmp2 $baseDir/align/$baseName.tmp3 > $baseDir/align/$baseName
    rm -f $baseDir/align/$baseName.tmp $baseDir/align/$baseName.tmp2 $baseDir/align/$baseName.tmp3
done

# generate reference file
echo "generate reference file"

awk 'BEGIN{OFS="\t";FS="\t";}
{print "M"$4, "microT"}' $baseDir/align/mm10_microT_sb3_align.bed > $baseDir/reference/mm10_bindingID_reference.txt
awk 'BEGIN{OFS="\t";FS="\t";}
{print "M"$4, "miRanda"}' $baseDir/align/mm10_miranda_cons_sb3_align.bed >> $baseDir/reference/mm10_bindingID_reference.txt
awk 'BEGIN{OFS="\t";FS="\t";}
{print "M"$4, "miRmap"}' $baseDir/align/mm10_miRmap_cons_sb3_align.bed >> $baseDir/reference/mm10_bindingID_reference.txt
awk 'BEGIN{OFS="\t";FS="\t";}
{print "M"$4, "PicTar"}' $baseDir/align/mm10_pictar_cons_sb3_align.bed >> $baseDir/reference/mm10_bindingID_reference.txt
awk 'BEGIN{OFS="\t";FS="\t";}
{print "M"$4, "PITA"}' $baseDir/align/mm10_PITA_cons_sb3_align.bed >> $baseDir/reference/mm10_bindingID_reference.txt
awk 'BEGIN{OFS="\t";FS="\t";}
{print "M"$4, "RNA22"}' $baseDir/align/mm10_RNA22_cons_sb3_align.bed >> $baseDir/reference/mm10_bindingID_reference.txt
awk 'BEGIN{OFS="\t";FS="\t";}
{print "M"$4, "TargetScan"}' $baseDir/align/mm10_targetscan_cons_sb3_align.bed >> $baseDir/reference/mm10_bindingID_reference.txt

# cat and sort predict files
# sorted by: miRNAID,chr,strand,start,end,siteID
echo "cat and sort predict files"
echo "sorted by: miRNAID,chr,strand,start,end,siteID"
cat $baseDir/align/mm10_microT_sb3_align.bed \
  $baseDir/align/mm10_miranda_cons_sb3_align.bed $baseDir/align/mm10_miRmap_cons_sb3_align.bed \
  $baseDir/align/mm10_pictar_cons_sb3_align.bed $baseDir/align/mm10_PITA_cons_sb3_align.bed \
  $baseDir/align/mm10_RNA22_cons_sb3_align.bed $baseDir/align/mm10_targetscan_cons_sb3_align.bed \
  | sort -t $'\t' -k7,7 -k1,1 -k6,6 -k2,2n -k3,3n -k4,4 | \
  awk 'BEGIN{OFS="\t";FS="\t";}{$4="M"$4;print $0;}' > $baseDir/align/mm10_all_sorted_cons_sb3_align.bed

endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s);
end_seconds=$(date --date="$endtime" +%s);
runtime=$((end_seconds-start_seconds));
echo "Total runtime：${runtime}s"

