#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
from functools import reduce
import datetime
from multiBioPro import BedMan

parser = argparse.ArgumentParser()
parser.add_argument('-anno', action='store', type=str,
                    help='The annotation bed12')
parser.add_argument('-input', action='store', type=str,
                    help='The input of original rriScan result file')
parser.add_argument('-prefix', action='store', type=str,
                    help='The prefix of output files')
parser.add_argument('-ref', action='store', type=str,
                    help='The reference files')
parser.add_argument('-output', action='store', type=str,
                    default="./", help='The output directory')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

def ReconPairing (leftLocus, rightLocus, leftPair, alignment, rightPair):
    leftPair = leftPair.replace("RNA1  5'>", '').replace(">3'", '')
    alignment = alignment.replace("Pair     ", '').replace('.', ' ')
    rightPair = rightPair.replace("RNA2  3'<", '').replace("<5'", '')
    ## remove the flanking unpairing base
    pairLen = len(alignment)
    trimStart = 0
    trimEnd = pairLen
    for i in range(pairLen):
        if alignment[i] != ' ':
            trimStart = i
            break
    for i in range(pairLen):
        if alignment[::-1][i] != ' ':
            trimEnd = i 
            break
    ## reconstruct pairing
    startIndex = trimStart
    endIndex = pairLen - trimEnd
    leftPairTrim = leftPair[startIndex:endIndex]
    alignmentTrim = alignment[startIndex:endIndex]
    rightPairTrim = rightPair[startIndex:endIndex]
    ## reconstruct pairing locus
    ## for gap in leftPair
    trimSeq = leftPair[0:trimStart].replace('-', '')
    trimSeqStart = len(trimSeq)
    trimSeq = leftPair[::-1][0:trimEnd].replace('-', '')
    trimSeqEnd = len(trimSeq)
    if leftLocus[3] == '+':
        leftLocus[1] = leftLocus[1] + trimSeqStart
        leftLocus[2] = leftLocus[2] - trimSeqEnd
    else:
        leftLocus[1] = leftLocus[1] + trimSeqEnd
        leftLocus[2] = leftLocus[2] - trimSeqStart
    ## for gap in rightPair
    trimSeq = rightPair[0:trimStart].replace('-', '')
    trimSeqStart = len(trimSeq)
    trimSeq = rightPair[::-1][0:trimEnd].replace('-', '')
    trimSeqEnd = len(trimSeq)
    if rightLocus[3] == '+':
        rightLocus[1] = rightLocus[1] + trimSeqEnd
        rightLocus[2] = rightLocus[2] - trimSeqStart
    else:
        rightLocus[1] = rightLocus[1] + trimSeqStart
        rightLocus[2] = rightLocus[2] - trimSeqEnd
    leftLocus = ':'.join(map(str, leftLocus))
    rightLocus = ':'.join(map(str, rightLocus))
    ## return results
    return [leftPairTrim, alignmentTrim, rightPairTrim, leftLocus, rightLocus]

def SplitGeneInfo (txInfoList, geneType, txGeneDict, geneDict):
    global regionDict
    txid = txInfoList[0].split('.')[0]
    if txid in txGeneDict:
        geneID = txGeneDict[txid]
    else:
        if re.match(r'^ENS', txid):
            ##### ENST00000497648.1|RABGGTB-214|ENSG00000137955.16|RABGGTB|protein_coding|exon-2
            geneID = txInfoList[2].split('.')[0]
            if geneID not in geneDict:
                return False
        else:
            geneID = txInfoList[2]
            geneName = txInfoList[3]
            if geneType == 'rRNA':
                #### rmsk995909|LSU-rRNA_Hsa|LSU-rRNA_Hsa-rRNA|LSU-rRNA_Hsa-rRNA|rRNA|exon-1
                #### NR_146151_blat|RNA45SN3_blat|NR_146151_blat|RNA45SN3_blat|rRNA|exon-1
                if bool(re.search(r'rmsk', txid)) is True:
                    geneID = txInfoList[2].split('_')[0]
                    geneName = txInfoList[3].split('_')[0]
                else:
                    geneName = txInfoList[3].replace('_blat', '')
            if geneType == 'CDBox':
                geneType = 'snoRNA'
            geneDict[geneID] = [geneName, geneType]
    ## get region type
    region = geneType in regionDict and regionDict[geneType] or 'Exon'
    return [geneID, region, geneDict]

def MaxConMatchNumCount (leftSeq, alignment, rightSeq, watsonCrickDict):
    # maximum contiguous mathch Watson-Crick pair
    lreqList = list(leftSeq)
    pairList = list(alignment)
    rreqList = list(rightSeq)
    if len(lreqList) == len(rreqList) and len(rreqList) == len(pairList):
        for i in range(len(lreqList)):
            lread = lreqList[i]
            rread = rreqList[i]
            pair = pairList[i]
            if pair == '|' and rread in watsonCrickDict[lread]:
                pairList[i] = "="
    else:
        sys.stderr.write("alignment and sequences should be mathched")
        sys.exit()
    pair = ''.join(pairList)
    alignList = re.sub(r'(\.+)|(\|+)|(\s+)', '.', pair).split('.')
    maxConMatchNum = max(map(lambda x:x.count('='), alignList))
    return maxConMatchNum

def CreatWatsonCrickDict (watsonList, crickList):
    if len(watsonList) != len(crickList):
        sys.stderr.write("watson and crick list should be in the same length")
        sys.exit()
    pairDict = defaultdict(dict)
    for i in range(len(watsonList)):
        wlist = [watsonList[i].upper(), watsonList[i].lower()]
        clist = [crickList[i].upper(), crickList[i].lower()]
        for j in range(len(wlist)):
            for k in range(len(clist)):
                pairDict[wlist[j]][clist[k]] = 1
                pairDict[clist[k]][wlist[j]] = 1
    return pairDict

watsonList = ['A', 'U', 'C']
crickList = ['T', 'A', 'G']
watsonCrickDict = CreatWatsonCrickDict(watsonList, crickList)

txGeneDict = defaultdict(dict)
geneDict = defaultdict(dict)
with open (args.anno, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        inforList = row[3].split('|')
        txID = inforList[0].split('.')[0]
        geneID = inforList[3]
        geneName = inforList[4]
        geneType = inforList[5]
        txGeneDict[txID] = geneID
        geneDict[geneID] = [geneName, geneType]

refDict = defaultdict(dict)
with open(args.ref, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        sraAccesion = row[9]
        datasetID = row[0]
        seqType = row[2]
        refDict[sraAccesion] = datasetID
        refDict[datasetID] = seqType

trnaTypeDict = {'tRNA':1, 'Mt_tRNA':1}
regionDict = {'cds':'CDS', 'intron':'Intron', 'exon':'Exon', 'utr3':'3\' UTR', 'utr5':'5\' UTR', 'CDBox': 'CD Box'}
interactID = 1
interactDict = defaultdict(dict)
with open (args.input, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        ## hES_SRR3404927, MouseBrainTissue_SRR3500451.extendedFrags, skip invalid samples
        sraAccesion = re.findall(r'SRR\d+', row[35])[0]
        if sraAccesion not in refDict:
            continue
        ## skip tRNA and Mt_tRNA
        leftGeneType = row[30].replace('Pol3_', '')
        rightGeneType = row[33].replace('Pol3_', '')
        ## skip tRNA
        if leftGeneType in trnaTypeDict or rightGeneType in trnaTypeDict:
            continue
        ## skip rRNA-rRNA interactions
        if leftGeneType == 'rRNA' and rightGeneType == 'rRNA':
            continue
        ## rriScanGroup717-11_9869573-13_lr
        ## groupid-isoform_readid-readNum_leftread
        datasetID = refDict[sraAccesion]
        idInfoList = re.split(r'_|-', row[3])
        grouid = idInfoList[0].replace('rriScanGroup', '')
        isoformid = idInfoList[1]
        readid = idInfoList[2]
        readNum = idInfoList[3]
        uniqReadID = '-'.join([datasetID, readid, readNum])
        ## [chrom, start, end, strand]
        leftLocus = [row[0], int(row[1]), int(row[2]), row[5]]
        rightLocus = [row[6], int(row[7]), int(row[8]), row[11]]
        ## base pairing
        leftPair = row[19]
        alignment = row[20]
        rightPair = row[21]
        leftPair, alignment, rightPair, leftLocus, rightLocus = ReconPairing(leftLocus, rightLocus, leftPair, alignment, rightPair)
        maxConMatchNum = MaxConMatchNumCount(leftPair, alignment, rightPair, watsonCrickDict)
        ##
        freeEnergy = row[17]
        smithWaterman = row[23]
        ## gene info
        ## left
        leftTxInfoList = row[28].split('|')
        geneInfoRes = SplitGeneInfo(leftTxInfoList, leftGeneType, txGeneDict, geneDict)
        if bool(geneInfoRes) is False:
            continue
        leftGeneID, leftRegion, geneDict = geneInfoRes
        ### right
        rightTxInfoList = row[31].split('|')
        geneInfoRes = SplitGeneInfo(rightTxInfoList, rightGeneType, txGeneDict, geneDict)
        if bool(geneInfoRes) is False:
            continue
        rightGeneID, rightRegion, geneDict = geneInfoRes
        ## reconstruct locus
        interactLocus = '|'.join([leftLocus, rightLocus])
        revInteractLocus = '|'.join([rightLocus, leftLocus])
        if interactLocus not in interactDict:
            if revInteractLocus in interactDict:
                interactDict[revInteractLocus]['readID'].append(uniqReadID)
            else:
                interactDict[interactLocus] = defaultdict(list)
                interactDict[interactLocus]['interactID'] = interactID
                interactDict[interactLocus]['gene'] = [leftGeneID, rightGeneID]
                interactDict[interactLocus]['region'] = [leftRegion, rightRegion]
                interactDict[interactLocus]['freeEnergy'] = freeEnergy
                interactDict[interactLocus]['leftPair'] = leftPair
                interactDict[interactLocus]['alignment'] = alignment
                interactDict[interactLocus]['rightPair'] = rightPair
                interactDict[interactLocus]['smithWaterman'] = smithWaterman
                interactDict[interactLocus]['maxConMatchNum'] = maxConMatchNum
                interactDict[interactLocus]['readID'].append(uniqReadID)
                interactID += 1
        else:
            interactDict[interactLocus]['readID'].append(uniqReadID)

geneInteractDict = defaultdict(list)
for interactLocus in interactDict.keys():
    locusList = [interactLocus, '']
    leftLocus, rightLocus = interactLocus.split('|')
    leftGeneID, rightGeneID = interactDict[interactLocus]['gene']
    if leftGeneID == rightGeneID:
        continue
    geneKey = '|'.join([leftGeneID, rightGeneID])
    revGeneKey = '|'.join([rightGeneID, leftGeneID])
    if geneKey not in geneInteractDict and revGeneKey not in geneInteractDict:
        geneInteractDict[geneKey].append(locusList)
    else:
        if geneKey in geneInteractDict:
            geneInteractDict[geneKey].append(locusList)
        else:
            revInteractLocus = '|'.join(interactLocus.split('|')[::-1])
            locusList = [revInteractLocus, interactLocus]
            geneInteractDict[revGeneKey].append(locusList)

rnaRNApairFile = os.path.join(args.output, args.prefix + '_rnaRNApair.txt')
with open(rnaRNApairFile, 'w') as out:
    tempList = ['interactID', 'geneID', 'geneName', 'geneType', 'pairGeneID', 'pairGeneName',
        'pairGeneType', 'interactLocus','expNum', 'seqTypeNum', 'readsNum', 'freeEnergy',
        'smithWaterman', 'maxConMatchNum', 'leftRegion', 'rightRegion', 'leftPair', 'alignment', 'rightPair', 'readIDcat']
    out.write('\t'.join(tempList) + '\n')
    for interactLocus in sorted(interactDict.keys(), key=lambda x:interactDict[x]['interactID']):
        tempDict = interactDict[interactLocus]
        interactID = tempDict['interactID']
        readIDlist = sorted(set(tempDict['readID']))
        readIDcat = ','.join(readIDlist)
        expNum = len(set(map(lambda x:x.split('-')[0], readIDlist)))
        seqTypeNum = len(set(map(lambda x:refDict[x.split('-')[0]], readIDlist)))
        readsNum = sum(map(lambda x:int(x.split('-')[2]), readIDlist))
        geneList = tempDict['gene']
        leftGeneID = geneList[0]
        leftGeneName = geneDict[leftGeneID][0]
        leftGeneType = geneDict[leftGeneID][1]
        rightGeneID = geneList[1]
        rightGeneName = geneDict[rightGeneID][0]
        rightGeneType = geneDict[rightGeneID][1]
        leftRegion, rightRegion = tempDict['region']
        freeEnergy = tempDict['freeEnergy']
        smithWaterman = tempDict['smithWaterman']
        maxConMatchNum = tempDict['maxConMatchNum']
        leftPair = tempDict['leftPair']
        alignment = tempDict['alignment']
        rightPair = tempDict['rightPair']
        tempList = [str(interactID), leftGeneID, leftGeneName, leftGeneType, rightGeneID, rightGeneType,
            rightGeneName, interactLocus, str(expNum), str(seqTypeNum), str(readsNum), freeEnergy,
            smithWaterman, str(maxConMatchNum), leftRegion, rightRegion, leftPair, alignment, rightPair, readIDcat]
        out.write('\t'.join(tempList) + '\n')

rnaRNAfile = os.path.join(args.output, args.prefix + '_rnaRNA.txt')
with open(rnaRNAfile, 'w') as out:
    tempList = ['lineID', 'geneID', 'geneName', 'geneType', 'pairGeneID', 'pairGeneName', 'pairGeneType']
    tempList.extend(['interactionNum', 'interactIDcat', 'expNum', 'seqTypeNum', 'readsNum', 'minFreeEnergy', 'maxSmithWaterman', 'maxConMatchNum'])
    out.write('\t'.join(tempList) + '\n')
    lineID = 1
    for geneKey in sorted(geneInteractDict.keys()):
        geneKeyList = geneKey.split('|')
        leftGeneID = geneKeyList[0]
        leftGeneName = geneDict[leftGeneID][0]
        leftGeneType = geneDict[leftGeneID][1]
        rightGeneID = geneKeyList[1]
        rightGeneName = geneDict[rightGeneID][0]
        rightGeneType = geneDict[rightGeneID][1]
        readIDtotalSet = set()
        freeEnergyList = list()
        smithWatermanList = list()
        maxConMatchNumList = list()
        interactIDlist = list()
        interactLocusList = geneInteractDict[geneKey]
        for locusList in interactLocusList:
            if locusList[1]:
                tempDict = interactDict[locusList[1]]
            else:
                tempDict = interactDict[locusList[0]]
            interactID = str(tempDict['interactID'])
            interactIDlist.append(interactID)
            readIDlist = tempDict['readID']
            readIDtotalSet.update(readIDlist)
            freeEnergyList.append(float(tempDict['freeEnergy']))
            smithWatermanList.append(float(tempDict['smithWaterman']))
            maxConMatchNumList.append(tempDict['maxConMatchNum'])

        interactIDcat = ','.join(sorted(set(interactIDlist)))
        expNum = len(set(map(lambda x:x.split('-')[0], readIDtotalSet)))
        seqTypeNum = len(set(map(lambda x:refDict[x.split('-')[0]], readIDtotalSet)))
        readsNum = sum(map(lambda x:int(x.split('-')[2]), readIDtotalSet))
        minFreeEnergy = min(freeEnergyList)
        maxSmithWaterman = max(smithWatermanList)
        maxConMatchNum = max(maxConMatchNumList)
        tempList = [str(lineID), leftGeneID, leftGeneName, leftGeneType, rightGeneID, rightGeneName, rightGeneType]
        tempList.extend([len(interactLocusList), interactIDcat, expNum, seqTypeNum, readsNum, minFreeEnergy, maxSmithWaterman, maxConMatchNum])
        out.write('\t'.join(list(map(str,tempList))) + '\n')
        lineID += 1
