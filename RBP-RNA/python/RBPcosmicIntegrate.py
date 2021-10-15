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
        geneDict[geneID] = geneName
        if geneID not in RBPgeneDict[RBP]:
            RBPgeneDict[RBP][geneID] = defaultdict(dict)
            RBPgeneDict[RBP][geneID][cosmicID] = defaultdict(set)
            RBPgeneDict[RBP][geneID][cosmicID]['clipID'].update([clipID])
            RBPgeneDict[RBP][geneID][cosmicID]['sample'].update(sampleList)
        else:
            if cosmicID not in RBPgeneDict[RBP][geneID]:
                RBPgeneDict[RBP][geneID][cosmicID] = defaultdict(set)
                RBPgeneDict[RBP][geneID][cosmicID]['clipID'].update([clipID])
                RBPgeneDict[RBP][geneID][cosmicID]['sample'].update(sampleList)
            else:
                RBPgeneDict[RBP][geneID][cosmicID]['clipID'].update([clipID])
                RBPgeneDict[RBP][geneID][cosmicID]['sample'].update(sampleList)

headerList = ['lineID', 'RBP', 'geneID', 'geneName', 'cosmicID',
    'tissueNum', 'diseaseNum', 'sampleNum', 'clipExpNum',
    'clipIDnum', 'sampleCat', 'clipIDcat']
RBPlist = sorted(RBPgeneDict.keys())
with open(args.output, 'w') as out:
    lineID = 1
    out.write('\t'.join(headerList) + '\n')
    for RBP in RBPlist:
        geneIDlist = sorted(RBPgeneDict[RBP].keys())
        for geneID in geneIDlist:
            geneName = geneDict[geneID]
            cosmicIDlist = sorted(RBPgeneDict[RBP][geneID].keys())
            for cosmicID in cosmicIDlist:
                tempDict = RBPgeneDict[RBP][geneID][cosmicID]
                sampleList = sorted(tempDict['sample'])
                sampleIDlist = list()
                tissueIDlist = list()
                diseaseIDlist = list()
                for sample in sampleList:
                    sampleID = sample.split('|')[0]
                    tissueID, diseaseID = sample.split('|')[1].split('-')
                    sampleIDlist.append(sampleID)
                    tissueIDlist.append(tissueID)
                    diseaseIDlist.append(diseaseID)
                clipIDlist = sorted(tempDict['clipID'])
                clipExpList = list()
                for clipID in clipIDlist:
                    expID = clipID.split('-')[0]
                    clipExpList.append(expID)
                sampleNum = str(len(set(sampleIDlist)))
                tissueNum = str(len(set(tissueIDlist)))
                diseaseNum = str(len(set(diseaseIDlist)))
                clipExpNum = str(len(set(clipExpList)))
                clipIDnum = str(len(set(clipIDlist)))
                sampleCat = ','.join(sampleList)
                clipIDcat = ','.join(clipIDlist)
                tempList = [str(lineID), RBP, geneID, geneName, cosmicID,
                    tissueNum, diseaseNum, sampleNum, clipExpNum,
                    clipIDnum, sampleCat, clipIDcat]
                out.write('\t'.join(tempList) + '\n')
                lineID += 1
