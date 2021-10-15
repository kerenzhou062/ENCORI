#!/usr/bin/env python3
import os
import sys
import argparse
import re
from collections import defaultdict
import datetime
from multiBioPro import BedMan

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The integrated mirTarget-seq file')
parser.add_argument('-output', action='store', type=str,
                    default="cluster.txt", help='The output cluster file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

starttime = datetime.datetime.now()

clusterDict = defaultdict(dict)
miRNAgeneExpDict = defaultdict(dict)
removeDict = defaultdict(list)
with open(args.input, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        cluster = row[3]
        narrowStart = int(row[1])
        narrowEnd = int(row[2])
        broadStart = int(row[8])
        broadEnd = int(row[9])
        bed12Row = row[40:]
        # [[exonblock], [intronblock], [thickup, thick, thickdown]]
        # or [[exonblock], [intronblock]]
        if narrowStart < int(bed12Row[1]) or narrowEnd > int(bed12Row[2]):
            removeDict['narrow'].append(cluster)
            continue
        if broadStart < int(bed12Row[1]) or broadEnd > int(bed12Row[2]):
            removeDict['broad'].append(cluster)
            continue
        decodeList = BedMan.decodeBed12(bed12Row)
        narrowLocus = [narrowStart, narrowEnd]
        if len(decodeList) == 2:
            continue
        txInfo = bed12Row[3].split('|')
        txID = txInfo[0]
        txName = txInfo[1]
        geneID = txInfo[3]
        geneName = txInfo[4]
        geneType = txInfo[5]
        strand = bed12Row[5]
        if strand == '+':
            thickTypeList = ['5\'UTR', 'CDS', '3\'UTR']
        else:
            thickTypeList = ['3\'UTR', 'CDS', '5\'UTR']
        exonBlock = decodeList[0]
        intronBlock = decodeList[1]
        thickBlock = decodeList[2]
        exonBlockOverlap = list(map(lambda x:BedMan.overlap(narrowLocus, x),
            exonBlock))
        thickBlockOverlapList = list()
        for i in thickBlock:
            if i:
                overlap = list(map(lambda x:BedMan.overlap(narrowLocus, x),
                    i))
                thickBlockOverlapList.append(overlap)
            else:
                thickBlockOverlapList.append([])
        exonOverlapList = list()
        exonCount = len(exonBlockOverlap)
        for i in range(exonCount):
            if exonBlockOverlap[i]:
                if strand == '+':
                    exonNum = 'Exon-' + str(i + 1)
                else:
                    exonNum = 'Exon-' + str(exonCount - i)
                exonOverlapList.append(exonNum)
        thickOverlapList = list()
        for i in range(3):
            if sum(thickBlockOverlapList[i]):
                thickOverlapList.append(thickTypeList[i])
        exonBlockLocate = ','.join(exonOverlapList)
        thickBlockLocate = ','.join(thickOverlapList)
        # cluster-gene statistics
        if cluster not in clusterDict:
            clusterDict[cluster]['info'] = row[0:40]
            clusterDict[cluster]['gene'] = defaultdict(list)
            clusterDict[cluster]['gene'][geneID].append([geneName, geneType, txID,
                txName, exonBlockLocate, thickBlockLocate])
        else:
            clusterDict[cluster]['gene'][geneID].append([geneName, geneType, txID,
                txName, exonBlockLocate, thickBlockLocate])
        # miRNA-gene statistics
        miRNAid = row[6]
        clipIDList = re.split(',|:', row[36])
        clipExpList = list(map(lambda x: x.split('-')[0], clipIDList))
        if row[39] != 'NA':
            degraIDList = clipDegraIDList = row[39].split(',')
            degraExpList = list(map(lambda x: x.split('-')[0], clipDegraIDList))
        else:
            degraExpList = list()
        if miRNAid not in miRNAgeneExpDict:
            miRNAgeneExpDict[miRNAid][geneID] = defaultdict(dict)
            miRNAgeneExpDict[miRNAid][geneID]['clip'] = clipExpList
            miRNAgeneExpDict[miRNAid][geneID]['degra'] = degraExpList
        else:
            if geneID not in miRNAgeneExpDict[miRNAid]:
                miRNAgeneExpDict[miRNAid][geneID] = defaultdict(dict)
                miRNAgeneExpDict[miRNAid][geneID]['clip'] = clipExpList
                miRNAgeneExpDict[miRNAid][geneID]['degra'] = degraExpList
            else:
                miRNAgeneExpDict[miRNAid][geneID]['clip'].extend(clipExpList)
                miRNAgeneExpDict[miRNAid][geneID]['degra'].extend(degraExpList)

for miRNAid in sorted(miRNAgeneExpDict.keys()):
    for geneID in sorted(miRNAgeneExpDict[miRNAid].keys()):
        miRNAgeneExpDict[miRNAid][geneID]['clip'] = str(len(set(miRNAgeneExpDict[miRNAid][geneID]['clip'])))
        miRNAgeneExpDict[miRNAid][geneID]['degra'] = str(len(set(miRNAgeneExpDict[miRNAid][geneID]['degra'])))

headerList = ['lineID', 'chromosome','narrowStart','narrowEnd','clusterID',
    'softwareCount','strand','miRNAid','miRNAname','broadStart','broadEnd', 'tdmdScore',
    'PITA', 'PicTar', 'RNA22', 'TargetScan', 'miRanda', 'miRmap', 'microT',
    'PITAID', 'PicTarID', 'RNA22ID', 'TargetScanID', 'miRandaID', 'miRmapID',
    'microTID', 'PITASeqStats', 'PicTarSeqStats', 'RNA22SeqStats',
    'TargetScanSeqStats', 'miRandaSeqStats', 'miRmapSeqStats', 'microTSeqStats',
    'clipExpNum','clusterClipExpNum','clipIDnum','RBP','RBPclipNum','RBPclipID',
    'degraExpNum', 'clusterDegraExpNum','degraIDNum','degraID', 'geneID', 'geneName',
    'geneType', 'txInfo']

with open(args.output, 'w') as out:
    out.write('\t'.join(headerList) + '\n')
    clusters = sorted(clusterDict.keys())
    count = 1
    for cluster in clusters:
        tempDict = clusterDict[cluster]
        infoList = tempDict['info']
        miRNAid = infoList[6]
        geneIDlist = sorted(tempDict['gene'].keys())
        for geneID in geneIDlist:
            # add clipExpNum and degraExpNum to clusterInfo
            tempInfoList = infoList[0:]
            clipExpNum = miRNAgeneExpDict[miRNAid][geneID]['clip']
            degraExpNum = miRNAgeneExpDict[miRNAid][geneID]['degra']
            tempInfoList[32] = '\t'.join([clipExpNum, tempInfoList[32]])
            tempInfoList[37] = '\t'.join([degraExpNum, tempInfoList[37]])
            tempInfo = '\t'.join(tempInfoList)
            tempTxList = tempDict['gene'][geneID]
            geneName = tempTxList[0][0]
            geneType = tempTxList[0][1]
            tempList = [geneID, geneName, geneType]
            # cat transcript info
            tempTxCatList = list()
            for tx in tempTxList:
                txInfo = ':'.join(tx[2:])
                tempTxCatList.append(txInfo)
            tempTx = '|'.join(tempTxCatList)
            tempList.append(tempTx)
            temp = '\t'.join(tempList) + '\n'
            out.write('\t'.join([str(count), tempInfo, temp]))
            count += 1

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
