#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
import datetime
from multiBioPro import BedMan

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The integrated mirTarget-seq file')
parser.add_argument('-cds', action='store_true',
                    default=False, help='cds flag')
parser.add_argument('-circRNA', action='store_true',
                    default=False, help='cds flag')
parser.add_argument('-type', action='store', type=str,
                    default='bed12', help='bed6 or bed12(default:bed12)')
parser.add_argument('-output', action='store', type=str,
                    default="RBPanno.txt", help='The output clusterID file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

RBPgeneDict = defaultdict(dict)
geneDict = defaultdict(dict)

with open(args.input, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        clusterID = row[4]
        RBP = row[9]
        clipIDList = row[12].split(',')
        geneID = row[13]
        geneName = row[14]
        geneType = row[15]
        if geneID not in RBPgeneDict[RBP]:
            RBPgeneDict[RBP][geneID] = defaultdict(set)
            RBPgeneDict[RBP][geneID]['clipID'].update(clipIDList)
            RBPgeneDict[RBP][geneID]['clusterID'].update([clusterID])
        else:
            RBPgeneDict[RBP][geneID]['clipID'].update(clipIDList)
            RBPgeneDict[RBP][geneID]['clusterID'].update([clusterID])
        geneDict[geneID] = [geneName, geneType]

headerList = ['lineID', 'RBP', 'geneID', 'geneName', 'geneType', 'clusterNum', 'clipExpNum', 'clipIDnum', 'clipID']
with open (args.output, 'w') as out:
    out.write('\t'.join(headerList) + '\n')
    lineID = 1
    for RBP in sorted(RBPgeneDict.keys()):
        for geneID in sorted(RBPgeneDict[RBP].keys()):
            clusterIDnum = str(len(RBPgeneDict[RBP][geneID]['clusterID']))
            clipIDlist = sorted(RBPgeneDict[RBP][geneID]['clipID'])
            clipExpIDlist = set(map(lambda x:x.split('-')[0], clipIDlist))
            clipExpNum = str(len(clipExpIDlist))
            clipIDnum = str(len(clipIDlist))
            clipIDcat = ','.join(clipIDlist)
            geneName = geneDict[geneID][0]
            geneType = geneDict[geneID][1]
            tempList = [str(lineID), RBP, geneID, geneName, geneType, clusterIDnum, clipExpNum, clipIDnum, clipIDcat]
            out.write('\t'.join(tempList) + '\n')
            lineID += 1
