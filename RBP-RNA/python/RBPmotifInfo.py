#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
import datetime
import math
import shutil
from glob import glob
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument('-input', action='store', type=str,
                    help='The input bed folder')
parser.add_argument('-meta', action='store', type=str,
                    help='metadata matrix of CLIP-seq')
parser.add_argument('-output', action='store', type=str,
                    default="cluster.txt", help='The output cluster file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

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
        seqType = row[2]
        cellTissue = row[4]
        metaDict[datasetID]['rbp'] = RBP
        metaDict[datasetID]['seqType'] = seqType
        metaDict[datasetID]['cellTissue'] = cellTissue


findComand = 'find "{0}" -type f -name "motif_information.txt"'.format(args.input)
findResult = bytes.decode(subprocess.check_output(findComand, shell=True))

infoFileList = list(map(os.path.realpath, findResult.rstrip().split('\n')))
motifDict = defaultdict(dict)
for infoFile in infoFileList:
    basenameList = re.split(r'\\|\/', infoFile)[1:]
    sampleID = basenameList[-2]
    motifDict[sampleID]['info'] = infoFile

findComand = 'find "{0}" -type f -name "motif_focus_information.txt"'.format(args.input)
findResult = bytes.decode(subprocess.check_output(findComand, shell=True))

focusFileList = list(map(os.path.realpath, findResult.rstrip().split('\n')))
for focusFile in focusFileList:
    basenameList = re.split(r'\\|\/', focusFile)[1:]
    sampleID = basenameList[-2]
    motifDict[sampleID]['focus'] = focusFile

lineID = 1
with open(args.output, 'w') as out:
    headerList = ['lineID', 'sampleID', 'rbp', 'seqType', 'cellTissue', 'motifID', 'motifSeq',
        'motifLength', 'totalPeakNum', 'totalBackgroundNum', 'pval', 'logPval',
        'targetPeakNum', 'percentTarget', 'backgroundPeakNum',
        'percentBackground', 'avPosition', 'bestMatch']
    out.write('\t'.join(headerList) + '\n')
    for sampleID in sorted(motifDict.keys()):
        rbp = metaDict[sampleID]['rbp']
        seqType = metaDict[sampleID]['seqType']
        cellTissue = metaDict[sampleID]['cellTissue']
        tempDict = defaultdict(dict)
        with open(motifDict[sampleID]['focus'], 'r') as f:
            for line in f:
                row = line.strip().split('\t')
                motifID = row[0]
                targetPeakNum = re.sub(r'\.\d+$', '', row[5])
                backgroundPeakNum = re.sub(r'\.\d+$', '', row[7])
                tempDict[motifID]['targetPeakNum'] = targetPeakNum
                tempDict[motifID]['backgroundPeakNum'] = backgroundPeakNum
        with open(motifDict[sampleID]['info'], 'r') as f:
            for line in f:
                row = line.strip().split('\t')
                motifID = row[0]
                motifSeq = row[1].replace('T', 'U')
                motifLength = str(len(motifSeq))
                totalPeakNum = row[2]
                totalBackgroundNum = row[3]
                pval = row[4]
                logPval = row[5]
                percentTarget = row[6].strip("%")
                percentBackground = row[7].strip("%")
                targetPeakNum = tempDict[motifID]['targetPeakNum']
                backgroundPeakNum = tempDict[motifID]['backgroundPeakNum']
                #Average Position of motif in Targets
                avPosition = row[8]
                bestMatch = row[9]
                tempList = [str(lineID), sampleID, rbp, seqType, cellTissue, motifID, motifSeq,
                    motifLength, totalPeakNum, totalBackgroundNum, pval, logPval,
                    targetPeakNum, percentTarget, backgroundPeakNum,
                    percentBackground, avPosition, bestMatch]
                out.write('\t'.join(tempList) + '\n')
                lineID += 1
