#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The input miR-anno file')
parser.add_argument('-type', action='store', type=str,
                    help='The input gene type (mRNA/lncRNA...)')
parser.add_argument('-pan', action='store', type=str,
                    help='The pancancer result file')
parser.add_argument('-output', action='store', type=str,
                    default="./annoPancancer.txt", help='The output anno-pancancer result')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

miRgenePanDict = defaultdict(dict)
with open(args.pan, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        miRNAid = row[0]
        coExpGeneID = row[2]
        pancancerCoExpList = row[4:]
        pancancerNum = 0
        for cancerCoExp in pancancerCoExpList:
            cancerResultList = cancerCoExp.split(',')
            coExpPearsonR = float(cancerResultList[0])
            coExpPval = float(cancerResultList[1])
            #coExpFDR = float(cancerResultList[2])
            if coExpPval <= 0.05 and coExpPearsonR < 0:
                pancancerNum += 1
        miRgenePanDict[miRNAid][coExpGeneID] = str(pancancerNum)

with open(args.input, 'r') as f, open(args.output, 'w') as out:
    headerList = f.readline().strip().split('\t')
    headerList.append('pancancerNum')
    out.write('\t'.join(headerList) + '\n')
    for line in f:
        row = line.strip().split('\t')
        miRNAid = row[7]
        if args.type == 'mRNA':
            geneID = row[43]
        else:
            geneID = row[23]
        try:
            pancancerNum = miRgenePanDict[miRNAid][geneID]
        except KeyError:
            pancancerNum = '0'
        row.append(pancancerNum)
        out.write('\t'.join(row) + '\n')

