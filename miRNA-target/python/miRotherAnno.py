#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
import datetime
from multiBioPro import BedMan

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The integrated mirTarget-seq file')
parser.add_argument('-circRNA', action='store_true',
                    default=False, help='circRNA flag')
parser.add_argument('-type', action='store', type=str,
                    default='bed12', help='bed6 or bed12(default:bed12)')
parser.add_argument('-output', action='store', type=str,
                    default="miRAnno.txt", help='The output siteID file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

starttime = datetime.datetime.now()

miRNAgeneExpDict = defaultdict(dict)
siteDict = defaultdict(dict)
if args.type == 'bed12':
    with open(args.input, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            siteID = row[3]
            start = int(row[1])
            end = int(row[2])
            bed12Row = row[20:]
            # [[exonblock], [intronblock], [thickup, thick, thickdown]]
            # or [[exonblock], [intronblock]]
            if start < int(bed12Row[1]) or end > int(bed12Row[2]):
                continue
            decodeList = BedMan.decodeBed12(bed12Row)
            locus = [start, end]
            if args.circRNA:
                pass
            else:
                if len(decodeList) == 3:
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
            exonBlockOverlap = list(map(lambda x:BedMan.overlap(locus, x),
                exonBlock))
            exonOverlapList = list()
            exonCount = len(exonBlockOverlap)
            for i in range(exonCount):
                if exonBlockOverlap[i]:
                    if strand == '+':
                        exonNum = 'Exon-' + str(i + 1)
                    else:
                        exonNum = 'Exon-' + str(exonCount - i)
                    exonOverlapList.append(exonNum)
            exonBlockLocate = ','.join(exonOverlapList)
            # design for circRNA
            if args.circRNA:
                if len(decodeList) == 3:
                    thickBlock = decodeList[2]
                    thickBlockOverlapList = list()
                    for i in thickBlock:
                        if i:
                            overlap = list(map(lambda x:BedMan.overlap(locus, x),i))
                            thickBlockOverlapList.append(overlap)
                        else:
                            thickBlockOverlapList.append([])
                    thickOverlapList = list()
                    for i in range(3):
                        if sum(thickBlockOverlapList[i]):
                            thickOverlapList.append(thickTypeList[i])
                    thickBlockLocate = ','.join(thickOverlapList)
                else:
                    thickBlockLocate = 'Exon'
            else:
                thickBlockLocate = 'Exon'
            if siteID not in siteDict:
                siteDict[siteID]['info'] = row[0:20]
                siteDict[siteID]['gene'] = defaultdict(list)
                siteDict[siteID]['gene'][geneID].append([geneName, geneType, txID,
                    txName, exonBlockLocate, thickBlockLocate])
            else:
                siteDict[siteID]['gene'][geneID].append([geneName, geneType, txID,
                    txName, exonBlockLocate, thickBlockLocate])
            # miRNA-gene statistics
            miRNAid = row[6]
            clipIDList = re.split(',|:', row[16])
            clipExpList = list(map(lambda x: x.split('-')[0], clipIDList))
            if row[19] != 'NA':
                degraIDList = clipDegraIDList = row[19].split(',')
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
else:
    with open(args.input, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            siteID = row[3]
            start = int(row[1])
            end = int(row[2])
            bed6Row = row[20:]
            if start < int(bed6Row[1]) or end > int(bed6Row[2]):
                continue
            geneInfo = bed6Row[3].split('|')
            geneID = geneInfo[3]
            geneName = geneInfo[4]
            geneType = geneInfo[5]
            if siteID not in siteDict:
                siteDict[siteID]['info'] = row[0:20]
                siteDict[siteID]['gene'] = list()
                siteDict[siteID]['gene'].append([geneID, geneName, geneType])
            else:
                siteDict[siteID]['gene'].append([geneID, geneName, geneType])
            # miRNA-gene statistics
            miRNAid = row[6]
            clipIDList = re.split(',|:', row[16])
            clipExpList = list(map(lambda x: x.split('-')[0], clipIDList))
            if row[19] != 'NA':
                degraIDList = clipDegraIDList = row[19].split(',')
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

if args.type == 'bed12':
    if args.circRNA:
        headerList = ['lineID', 'chromosome', 'start', 'end', 'bindID','merClass',
            'strand', 'miRNAid', 'miRNAname', 'miRseq', 'align', 'targetSeq', 'tdmdScore',
            'clipExpNum', 'bindClipExpNum', 'clipIDnum', 'RBP', 'RBPclipNum',
            'RBPclipID', 'degraExpNum', 'bindDegraExpNum', 'degraIDNum', 'degraID',
            'geneID', 'geneName', 'geneType', 'txIDcat', 'txInfo']
    else:
        headerList = ['lineID', 'chromosome', 'start', 'end', 'bindID','merClass',
            'strand', 'miRNAid', 'miRNAname', 'miRseq', 'align', 'targetSeq', 'tdmdScore',
            'clipExpNum', 'bindClipExpNum', 'clipIDnum', 'RBP', 'RBPclipNum',
            'RBPclipID', 'degraExpNum', 'bindDegraExpNum', 'degraIDNum', 'degraID',
            'geneID', 'geneName', 'geneType', 'txInfo']
else:
    headerList = ['lineID', 'chromosome', 'start', 'end', 'bindID','merClass',
        'strand', 'miRNAid', 'miRNAname', 'miRseq', 'align', 'targetSeq', 'tdmdScore',
        'clipExpNum', 'bindClipExpNum', 'clipIDnum', 'RBP', 'RBPclipNum',
        'RBPclipID', 'degraExpNum', 'bindDegraExpNum', 'degraIDNum', 'degraID',
        'geneID', 'geneName', 'geneType']

if args.type == 'bed12':
    with open(args.output, 'w') as out:
        out.write('\t'.join(headerList) + '\n')
        siteIDs = sorted(siteDict.keys())
        count = 1
        for siteID in siteIDs:
            tempDict = siteDict[siteID]
            infoList = tempDict['info']
            miRNAid = infoList[6]
            geneIDlist = sorted(tempDict['gene'].keys())
            for geneID in geneIDlist:
                # add clipExpNum and degraExpNum to clusterInfo
                tempInfoList = infoList[0:]
                clipExpNum = miRNAgeneExpDict[miRNAid][geneID]['clip']
                degraExpNum = miRNAgeneExpDict[miRNAid][geneID]['degra']
                tempInfoList[12] = '\t'.join([clipExpNum, tempInfoList[12]])
                tempInfoList[17] = '\t'.join([degraExpNum, tempInfoList[17]])
                tempInfo = '\t'.join(tempInfoList)
                # cat transcript info
                tempTxList = tempDict['gene'][geneID]
                geneName = tempTxList[0][0]
                geneType = tempTxList[0][1]
                tempList = [geneID, geneName, geneType]
                tempTxCatList = list()
                tempTxInfoCatList = list()
                for tx in tempTxList:
                    txInfo = ':'.join(tx[1:])
                    tempTxInfoCatList.append(txInfo)
                    txID = tx[2]
                    tempTxCatList.append(txID)
                if args.circRNA:
                    tempTxCat = ','.join(sorted(set(tempTxCatList)))
                    tempList.append(tempTxCat)
                tempTx = '|'.join(tempTxInfoCatList)
                tempList.append(tempTx)
                temp = '\t'.join(tempList) + '\n'
                out.write('\t'.join([str(count), tempInfo, temp]))
                count += 1
else:
    with open(args.output, 'w') as out:
        out.write('\t'.join(headerList) + '\n')
        siteIDs = sorted(siteDict.keys())
        count = 1
        for siteID in siteIDs:
            tempDict = siteDict[siteID]
            infoList = tempDict['info']
            miRNAid = infoList[6]
            for geneInfo in sorted(tempDict['gene'], key=lambda x:x[0]):
                geneID = geneInfo[0]
                geneName = geneInfo[1]
                geneType = geneInfo[2]
                # add clipExpNum and degraExpNum to clusterInfo
                tempInfoList = infoList[0:]
                clipExpNum = miRNAgeneExpDict[miRNAid][geneID]['clip']
                degraExpNum = miRNAgeneExpDict[miRNAid][geneID]['degra']
                tempInfoList[12] = '\t'.join([clipExpNum, tempInfoList[12]])
                tempInfoList[17] = '\t'.join([degraExpNum, tempInfoList[17]])
                tempInfo = '\t'.join(tempInfoList)
                out.write('\t'.join([str(count), tempInfo, geneID, geneName, geneType]) + '\n')
                count += 1

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
