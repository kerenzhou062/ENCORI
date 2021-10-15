#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-anno', action='store', nargs='+', type=str,
                    help='The annotation files')
parser.add_argument('-input', action='store', type=str,
                    help='The input rnaRNA-result file')
parser.add_argument('-output', action='store', type=str,
                    default="./rnaRNA_gene.txt", help='The output rnaRNA_gene result')

args = parser.parse_args()
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()

extAnnoDict = {}
extAnnoDict['protein_coding'] = 'mRNA'
extAnnoDict['lncRNA'] = 'lncRNA'
extAnnoDict['miRNA'] = 'miRNA'
extAnnoDict['misc_RNA'] = 'sncRNA'
extAnnoDict['Mt_rRNA'] = 'sncRNA'
extAnnoDict['Mt_tRNA'] = 'sncRNA'
extAnnoDict['4.5S_scRNA'] = 'sncRNA'
extAnnoDict['5S_rRNA'] = 'sncRNA'
extAnnoDict['5S_rRNA_pseudogene'] = 'sncRNA'
extAnnoDict['7SK_snRNA'] = 'sncRNA'
extAnnoDict['7SL_srpRNA'] = 'sncRNA'
extAnnoDict['BC200_scRNA'] = 'sncRNA'
extAnnoDict['RMRP_ribozyme'] = 'sncRNA'
extAnnoDict['RPPH1_ribozyme'] = 'sncRNA'
extAnnoDict['SNAR_snRNA'] = 'sncRNA'
extAnnoDict['tRNA'] = 'sncRNA'
extAnnoDict['U6_snRNA'] = 'sncRNA'
extAnnoDict['vaultRNA'] = 'sncRNA'
extAnnoDict['YRNA_scRNA'] = 'sncRNA'
extAnnoDict['ribozyme'] = 'sncRNA'
extAnnoDict['rRNA'] = 'sncRNA'
extAnnoDict['U6_snRNA'] = 'sncRNA'
extAnnoDict['scaRna'] = 'sncRNA'
extAnnoDict['scaRNA'] = 'sncRNA'
extAnnoDict['scRNA'] = 'sncRNA'
extAnnoDict['snoRNA'] = 'sncRNA'
extAnnoDict['snRNA'] = 'sncRNA'
extAnnoDict['pseudogene'] = 'pseudogene'
extAnnoDict['TR_pseudogene'] = 'pseudogene'
extAnnoDict['IG_pseudogene'] = 'pseudogene'
extAnnoDict['rRNA_pseudogene'] = 'pseudogene'
extAnnoDict['TEC'] = 'mRNA'
extAnnoDict['TR_gene'] = 'mRNA'
extAnnoDict['IG_gene'] = 'mRNA'

geneInfoDict = defaultdict(dict)
annoTypeList = ['mRNA', 'lncRNA', 'pseudogene', 'sncRNA']
for i in range(len(args.anno)):
    with open(args.anno[i], 'r') as f:
        for line in f:
            row = line.split('\t')
            info = row[3].split('|')
            geneID = info[3]
            geneInfoDict[geneID] = annoTypeList[i]

rnaRNAtypeAnnoDict = defaultdict(dict)
with open(args.input, 'r') as f:
    __ = f.readline()
    for line in f:
        row = line.strip().split("\t")
        geneID = row[1]
        geneName = row[2]
        geneType = row[3]
        pairGeneID = row[4]
        pairGeneName = row[5]
        pairGeneType = row[6]
        if geneID not in geneInfoDict:
            annoGeneType = extAnnoDict[geneType]
        else:
            annoGeneType = geneInfoDict[geneID]
        if pairGeneID not in geneInfoDict:
            annoPairGeneType = extAnnoDict[pairGeneType]
        else:
            annoPairGeneType = geneInfoDict[pairGeneID]
        key = ','.join([annoGeneType,geneID,geneName])
        rnaRNAtypeAnnoDict[key] = 1
        key = ','.join([annoPairGeneType,pairGeneID,pairGeneName])
        rnaRNAtypeAnnoDict[key] = 1

keyList = sorted(rnaRNAtypeAnnoDict.keys())
lineID = 1
with open(args.output, 'w') as out:
    out.write('lineID\tannoType\tgeneID\tgeneName\n');
    for key in keyList:
        tempList = key.split(',')
        annoType = tempList[0]
        geneID = tempList[1]
        geneName = tempList[2]
        out.write('\t'.join([str(lineID), annoType, geneID, geneName]) + '\n')
        lineID += 1

