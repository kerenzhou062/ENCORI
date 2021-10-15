#!/usr/bin/env python3
import os
import sys
import re
import argparse
import copy
from collections import defaultdict
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The input RBP-cosmic file')
parser.add_argument('-diseaseRef', action='store', type=str,
                    help='The reference file of cosmic disease')
parser.add_argument('-agoRef', action='store', type=str,
                    help='The AGO reference file')
parser.add_argument('-meta', action='store', type=str,
                    help='metadata matrix of CLIP-seq')
parser.add_argument('-output', action='store', type=str,
                    default="./hg19_cosmic_stats.txt", help='The output anno-pancancer result')
args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

agoRefDict = defaultdict(int)
if os.path.getsize(args.agoRef):
    with open(args.agoRef, 'r') as f:
        for line in f:
            row = line.strip()
            agoRefDict[row] = 1

metaDict = defaultdict(dict)
with open(args.meta, 'r') as f:
    # dataSetId,Species,SeqType,GeneSymbol,Cell/Tissue,Treatment,Source,
    # Accession-GSE,Accession,FileName-RBP,Assembly,newFileName,
    # CitationForShort,Citation,PubMedID,Title
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        datasetID = row[0]
        RBP = row[3]
        if RBP not in agoRefDict:
            metaDict[datasetID] = RBP

tissueList = list()
with open(args.diseaseRef, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        tissueID = row[1]
        tissueList.append(tissueID)
tissueList = sorted(set(tissueList), key=lambda x:int(x))
tissueNum = len(tissueList)

tissueDict = defaultdict(dict)
for tissueID in tissueList:
    tissueDict[tissueID] = defaultdict(dict)

geneDict = defaultdict(dict)
RBPgeneDict = defaultdict(dict)
with open(args.input, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        datasetID = row[6]
        if datasetID not in metaDict:
            continue
        RBP = metaDict[datasetID]
        clipID = row[3]
        cosmicID = row[10]
        geneID = row[13]
        geneName = row[14]
        sampleList = row[17].split(',')
        mutType = row[20]
        geneDict[geneID] = geneName
        if geneID not in RBPgeneDict[RBP]:
            RBPgeneDict[RBP][geneID] = defaultdict(dict)
            RBPgeneDict[RBP][geneID][cosmicID] = defaultdict(set)
            RBPgeneDict[RBP][geneID][cosmicID]['clipID'].update([clipID])
            RBPgeneDict[RBP][geneID][cosmicID]['sample'].update(sampleList)
            RBPgeneDict[RBP][geneID][cosmicID]['mutType'].update([mutType])
        else:
            if cosmicID not in RBPgeneDict[RBP][geneID]:
                RBPgeneDict[RBP][geneID][cosmicID] = defaultdict(set)
                RBPgeneDict[RBP][geneID][cosmicID]['clipID'].update([clipID])
                RBPgeneDict[RBP][geneID][cosmicID]['sample'].update(sampleList)
                RBPgeneDict[RBP][geneID][cosmicID]['mutType'].update([mutType])
            else:
                RBPgeneDict[RBP][geneID][cosmicID]['clipID'].update([clipID])
                RBPgeneDict[RBP][geneID][cosmicID]['sample'].update(sampleList)
                RBPgeneDict[RBP][geneID][cosmicID]['mutType'].update([mutType])

headerList = ['lineID', 'RBP', 'geneID', 'geneName', 'tissueIDcat', 'mutTypeIDcat',
    'totalTssueNum', 'totalMutTypeNum']
suffixList = ['diseaseNum', 'mutTypeNum', 'cosmicNum', 'sampleNum', 'clipExpNum','clipNum']
for suffix in suffixList:
    for tissueID in tissueList:
        headerList.append(tissueID + '_' + suffix)
for tissueID in tissueList:
    tempHeadList = ['diseaseCat', 'disSampleNumCat', 'disClipExpNumCat',
        'disClipNumCat', 'disCosCat', 'disCosSampleNumCat',
        'cosmicCat', 'cosSampleNumCat', 'cosClipExpNumCat',
        'cosClipNumCat', 'cosDisCat', 'cosDisSampleNumCat']
    tempList = list(map(lambda x: tissueID + '_' + x, tempHeadList))
    headerList.extend(tempList)

RBPlist = sorted(RBPgeneDict.keys())
with open(args.output, 'w') as out:
    lineID = 1
    out.write('\t'.join(headerList) + '\n')
    for RBP in RBPlist:
        geneIDlist = sorted(RBPgeneDict[RBP].keys())
        for geneID in geneIDlist:
            geneName = geneDict[geneID]
            totalTissueIDlist = list()
            totalMutTypeIDlist = list()
            diseaseNumList = list()
            cosmicNumList = list()
            mutTypeNumList = list()
            sampleNumList = list()
            clipExpNumList = list()
            clipNumList = list()
            tissueInfo = list()
            diseaseDict = copy.deepcopy(tissueDict)
            cosmicDict = copy.deepcopy(tissueDict)
            cosmicIDlist = sorted(RBPgeneDict[RBP][geneID].keys())
            for cosmicID in cosmicIDlist:
                tempDict = RBPgeneDict[RBP][geneID][cosmicID]
                sampleList = sorted(tempDict['sample'])
                clipIDlist = sorted(tempDict['clipID'])
                mutTypeIDlist = sorted(tempDict['mutType'])
                totalMutTypeIDlist.extend(mutTypeIDlist)
                for sample in sampleList:
                    sampleID = sample.split('|')[0]
                    tissueID, diseaseID = sample.split('|')[1].split('-')
                    if diseaseID not in diseaseDict:
                        diseaseDict[tissueID][diseaseID][cosmicID] = defaultdict(set)
                        diseaseDict[tissueID][diseaseID][cosmicID]['sample'].update([sampleID])
                        diseaseDict[tissueID][diseaseID][cosmicID]['clipID'].update(clipIDlist)
                        diseaseDict[tissueID][diseaseID][cosmicID]['mutType'].update(mutTypeIDlist)
                    else:
                        if cosmicID not in diseaseDict[tissueID][diseaseID]:
                            diseaseDict[tissueID][diseaseID][cosmicID] = defaultdict(set)
                            diseaseDict[tissueID][diseaseID][cosmicID]['sample'].update([sampleID])
                            diseaseDict[tissueID][diseaseID][cosmicID]['clipID'].update(clipIDlist)
                            diseaseDict[tissueID][diseaseID][cosmicID]['mutType'].update(mutTypeIDlist)
                        else:
                            diseaseDict[tissueID][diseaseID][cosmicID]['sample'].update([sampleID])
                            diseaseDict[tissueID][diseaseID][cosmicID]['clipID'].update(clipIDlist)
                            diseaseDict[tissueID][diseaseID][cosmicID]['mutType'].update(mutTypeIDlist)
                    if cosmicID not in cosmicDict:
                        cosmicDict[tissueID][cosmicID][diseaseID] = defaultdict(set)
                        cosmicDict[tissueID][cosmicID][diseaseID]['sample'].update([sampleID])
                        cosmicDict[tissueID][cosmicID][diseaseID]['clipID'].update(clipIDlist)
                    else:
                        if diseaseID not in cosmicDict[tissueID][cosmicID]:
                            cosmicDict[tissueID][cosmicID][diseaseID] = defaultdict(set)
                            cosmicDict[tissueID][cosmicID][diseaseID]['sample'].update([sampleID])
                            cosmicDict[tissueID][cosmicID][diseaseID]['clipID'].update(clipIDlist)
                        else:
                            cosmicDict[tissueID][cosmicID][diseaseID]['sample'].update([sampleID])
                            cosmicDict[tissueID][cosmicID][diseaseID]['clipID'].update(clipIDlist)
            for tissueID in tissueList:
                if bool(diseaseDict[tissueID]):
                    totalTissueIDlist.append(tissueID)
                    mutTypeList = list()
                    cosmicList = list()
                    sampleList = list()
                    clipList = list()
                    # category by disease
                    diseaseList = sorted(diseaseDict[tissueID].keys())
                    disCosList = list()
                    disSampleNumList = list()
                    disClipExpNumList = list()
                    disClipNumList = list()
                    disCosSampleNumCatList = list()
                    for diseaseID in diseaseList:
                        disSampleList = list()
                        disClipList = list()
                        cosmicIDList = sorted(diseaseDict[tissueID][diseaseID].keys())
                        disCosSampleNumList = list()
                        disCosClipNumList = list()
                        for cosmicID in cosmicIDList:
                            tempDict = diseaseDict[tissueID][diseaseID][cosmicID]
                            mutTypeList.extend(tempDict['mutType'])
                            sampleList.extend(tempDict['sample'])
                            clipList.extend(tempDict['clipID'])
                            disSampleList.extend(tempDict['sample'])
                            disClipList.extend(tempDict['clipID'])
                            disCosSampleNum = str(len(tempDict['sample']))
                            disCosClipNum = str(len(tempDict['clipID']))
                            disCosSampleNumList.append(disCosSampleNum)
                            disCosClipNumList.append(disCosClipNum)
                        cosmicCat = '|'.join(cosmicIDList)
                        disCosSampleNumCat = '|'.join(disCosSampleNumList)
                        disSampleNumList.append(str(len(set(disSampleList))))
                        disClipExpNumList.append(str(len(set(map(
                            lambda x:x.split('-')[0], disClipList)))))
                        disClipNumList.append(str(len(set(disClipList))))
                        disCosList.append(cosmicCat)
                        disCosSampleNumCatList.append(disCosSampleNumCat)
                    diseaseCat = ','.join(diseaseList)
                    disSampleNumCat = ','.join(disSampleNumList)
                    disClipExpNumCat = ','.join(disClipExpNumList)
                    disClipNumCat = ','.join(disClipNumList)
                    disCosCat = ','.join(disCosList)
                    disCosSampleNumCat = ','.join(disCosSampleNumCatList)
                    tissueInfo.extend([diseaseCat, disSampleNumCat, disClipExpNumCat, disClipNumCat, disCosCat, disCosSampleNumCat])
                    # category by cosmicID
                    cosmicList = sorted(cosmicDict[tissueID].keys())
                    cosDisList = list()
                    cosSampleNumList = list()
                    cosClipExpNumList = list()
                    cosClipNumList = list()
                    cosDisSampleNumCatList = list()
                    for cosmicID in cosmicList:
                        cosSampleList = list()
                        cosClipList = list()
                        dieaseIDList = sorted(cosmicDict[tissueID][cosmicID].keys())
                        cosDisSampleNumList = list()
                        for diseaseID in dieaseIDList:
                            tempDict = cosmicDict[tissueID][cosmicID][diseaseID]
                            cosSampleList.extend(tempDict['sample'])
                            cosClipList.extend(tempDict['clipID'])
                            cosDisSampleNum = str(len(tempDict['sample']))
                            cosDisSampleNumList.append(cosDisSampleNum)
                        diseaseCat = '|'.join(dieaseIDList)
                        cosDisSampleNumCat = '|'.join(cosDisSampleNumList)
                        cosSampleNumList.append(str(len(set(cosSampleList))))
                        cosClipExpNumList.append(str(len(set(map(
                            lambda x:x.split('-')[0], cosClipList)))))
                        cosClipNumList.append(str(len(set(cosClipList))))
                        cosDisList.append(diseaseCat)
                        cosDisSampleNumCatList.append(cosDisSampleNumCat)
                    cosmicCat = ','.join(cosmicList)
                    cosSampleNumCat = ','.join(cosSampleNumList)
                    cosClipExpNumCat = ','.join(cosClipExpNumList)
                    cosClipNumCat = ','.join(cosClipNumList)
                    cosDisCat = ','.join(cosDisList)
                    cosDisSampleNumCat = ','.join(cosDisSampleNumCatList)
                    tissueInfo.extend([cosmicCat, cosSampleNumCat, cosClipExpNumCat, cosClipNumCat, cosDisCat, cosDisSampleNumCat])
                    diseaseNum = str(len(set(diseaseList)))
                    mutTypeNum = str(len(set(mutTypeList)))
                    cosmicNum = str(len(set(cosmicList)))
                    sampleNum = str(len(set(sampleList)))
                    clipExpNum = str(len(set(map(
                        lambda x:x.split('-')[0], clipList))))
                    clipNum = str(len(set(clipList)))
                else:
                    diseaseNum = '0'
                    mutTypeNum = '0'
                    cosmicNum = '0'
                    sampleNum = '0'
                    clipExpNum = '0'
                    clipNum = '0'
                    tissueInfo.extend(['N' for i in range(12)])
                diseaseNumList.append(diseaseNum)
                mutTypeNumList.append(mutTypeNum)
                cosmicNumList.append(cosmicNum)
                sampleNumList.append(sampleNum)
                clipExpNumList.append(clipExpNum)
                clipNumList.append(clipNum)
            tissueNum = len(totalTissueIDlist)
            totalTissueIDcat = ','.join(totalTissueIDlist)
            totalMutTypeNum = len(set(totalMutTypeIDlist))
            totalMutTypeIDcat = ','.join(sorted(set(totalMutTypeIDlist)))
            tempList = [str(lineID), RBP, geneID, geneName, totalTissueIDcat, totalMutTypeIDcat, str(tissueNum), str(totalMutTypeNum)]
            tempList.extend(diseaseNumList)
            tempList.extend(mutTypeNumList)
            tempList.extend(cosmicNumList)
            tempList.extend(sampleNumList)
            tempList.extend(clipExpNumList)
            tempList.extend(clipNumList)
            tempList.extend(tissueInfo)
            out.write('\t'.join(tempList) + '\n')
            lineID += 1
