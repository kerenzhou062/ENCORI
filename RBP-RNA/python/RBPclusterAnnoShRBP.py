#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
import datetime
from glob import glob
import math

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The integrated mirTarget-seq file')
parser.add_argument('-shFolder', action='store', type=str,
                    help='The folder containing shRBP data')
parser.add_argument('-ref', action='store', type=str,
                    help='The RBP reference')
parser.add_argument('-stats', action='store_true',
                    help='The stats mode')
parser.add_argument('-output', action='store', type=str,
                    default="RBPanno.txt", help='The output clusterID file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

RBPrefDict = defaultdict(str)
with open(args.ref, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        #RBPname,geneID,RBPOfficial
        RBPrefDict[row[0]] = row[2]

shFdrFiles = glob(os.path.join(args.shFolder, '*.ebseqFdr'))

cellLineList = list()
shRBPdict = defaultdict(dict)
shRBPexpDict = defaultdict(set)
for shFdrFile in shFdrFiles:
    basename = os.path.splitext(os.path.basename(shFdrFile))[0]
    nameContentList = basename.split('_')
    cellLine = nameContentList[0]
    cellLineList.append(cellLine)
    RBP = nameContentList[1]
    expID = nameContentList[2]
    controlID = nameContentList[3]
    shRBPexpDict[RBP].update([expID, controlID])
    with open(shFdrFile, 'r') as f:
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            geneID = row[0].split('.')[0].replace('"', '')
            postFClog2 = ("%.3f" % math.log(float(row[3]),2))
            if geneID not in shRBPdict[RBP]:
                shRBPdict[RBP][geneID] = defaultdict(str)
                shRBPdict[RBP][geneID][cellLine] = postFClog2
            else:
                shRBPdict[RBP][geneID][cellLine] = postFClog2

statsDict = defaultdict(dict)
expSet = set()
cellLineSet = sorted(set(cellLineList))
if args.stats:
    with open(args.input, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            RBP = row[1]
            RBPname = RBPrefDict[RBP]
            geneID = row[2]
            if RBPname in shRBPdict:
                expSet.update(list(shRBPexpDict[RBPname]))
                if geneID in shRBPdict[RBPname]:
                    tempDict = shRBPdict[RBPname][geneID]
                    for cellLine in cellLineSet:
                        if cellLine in tempDict:
                            statsDict['RBP-cell'][RBPname+'-'+cellLine] = 1
    print("shRNA screen data (RBP-cell line) number:", len(statsDict['RBP-cell'].keys()))
    print("shRNA screen data (experiment) number:", len(expSet) * 2)
else:
    with open(args.input, 'r') as f, open(args.output, 'w') as out:
        headerList = f.readline().strip().split('\t')
        headerList.extend(cellLineSet)
        out.write('\t'.join(headerList) + '\n')
        for line in f:
            row = line.strip().split('\t')
            RBP = row[1]
            RBPname = RBPrefDict[RBP]
            geneID = row[2]
            tempDiffList = list()
            if RBPname in shRBPdict:
                if geneID in shRBPdict[RBPname]:
                    tempDict = shRBPdict[RBPname][geneID]
                    for cellLine in cellLineSet:
                        if cellLine in tempDict:
                            tempDiffList.append(tempDict[cellLine])
                        else:
                            tempDiffList.append('0')
                else:
                    tempDiffList = ['0', '0']
            else:
                tempDiffList = ['0', '0']
            row.extend(tempDiffList)
            out.write('\t'.join(row) + '\n')
