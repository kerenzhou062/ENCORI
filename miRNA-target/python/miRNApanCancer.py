#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
import datetime
from glob import glob
import tempfile
import subprocess
from multiprocessing import Pool, Manager
import math

parser = argparse.ArgumentParser()
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-cutoff', action='store', type=float,
                    default=1, help='The miRNA reference')
parser.add_argument('-input', action='store', type=str,
                    help='The input miRNA-geneID or miRNA-geneID file')
parser.add_argument('-miRNAseq', action='store', type=str,
                    help='The folder contain all pancancer miRNA-seq files')
parser.add_argument('-RNAseq', action='store', type=str,
                    help='The folder contain all pancancer RNA-seq files')
parser.add_argument('-output', action='store', type=str,
                    default="./pancancer.txt", help='The output pancancer result')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def runCoExp(miRNAid, geneList):
    coExpDict = defaultdict(dict)
    coExpGeneDict = defaultdict(dict)
    for cancer in cancerList:
        miRNAseqFile = cancerDict[cancer]['miRNA']
        RNAseqFile = cancerDict[cancer]['RNA']
        geneExpDict = defaultdict(dict)
        with open(miRNAseqFile, 'r') as f:
            __ = f.readline()
            for line in f:
                row = line.strip().split('\t')
                geneID = row[0]
                if geneID == miRNAid:
                    miRNAname = miRNAnameDict[miRNAid]
                    geneExpTempList = row[2].split(',')
                    miRNAExplist = list(map(lambda x:math.log(float(x) + args.cutoff, 2), geneExpTempList))
        with open(RNAseqFile, 'r') as f:
            __ = f.readline()
            for line in f:
                row = line.strip().split('\t')
                geneID = row[0]
                if geneID in geneList:
                    geneExpDict[geneID]['name'] = row[1]
                    geneExpTempList = row[2].split(',')
                    geneExpDict[geneID]['exp'] = list(map(lambda x:math.log(float(x) + args.cutoff, 2), geneExpTempList))

        with tempfile.NamedTemporaryFile() as temp:
            miRNACoExpRow = [miRNAid, miRNAname]
            miRNACoExpRow.extend(miRNAExplist)
            tempLine = '\t'.join(list(map(str, miRNACoExpRow))) + '\n'
            temp.write(tempLine.encode('utf-8'))
            for geneID in geneList:
                if geneID in geneExpDict:
                    geneCoExpRow = [geneID, geneExpDict[geneID]['name']]
                    geneCoExpRow.extend(geneExpDict[geneID]['exp'])
                    tempLine = '\t'.join(list(map(str, geneCoExpRow))) + '\n'
                    temp.write(tempLine.encode('utf-8'))
            temp.flush()
            coExpCommand = 'coExpressionFDR -n 1 -p 10 -q 10 {0}'.format(temp.name)
            coExpResult = bytes.decode(subprocess.check_output(coExpCommand, shell=True))
            coExpResultList = coExpResult.split('\n')[1:]
            for line in coExpResultList:
                if line == '':
                    continue
                row = line.strip().split('\t')
                coExpGeneID = row[3]
                coExpGeneName = row[4]
                coExpPearsonR = row[6]
                coExpPval = row[7]
                coExpFDR = row[11]
                coExpGeneDict[coExpGeneID] = coExpGeneName
                coExpDict[coExpGeneID][cancer] = [coExpPearsonR, coExpPval, coExpFDR]
    return [miRNAid, coExpGeneDict, coExpDict]

starttime = datetime.datetime.now()
sys.stderr.write("Program starting!\n")

miRNAnameDict = defaultdict(str)
pairDict = defaultdict(list)
with open(args.input, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        miRNAnameDict[row[0]] = row[1]
        if row[2] not in pairDict:
            pairDict[row[0]].append(row[2])

cancerDict = defaultdict(dict)
miRNAseqFiles = glob(os.path.join(args.miRNAseq, '*_miRNAseq.txt'))

for miRNAseqFile in miRNAseqFiles:
    cancerName = os.path.splitext(os.path.basename(miRNAseqFile))[0].split('_')[0]
    cancerDict[cancerName]['miRNA'] = miRNAseqFile

RNAseqFiles = glob(os.path.join(args.RNAseq, '*_RNAseq.txt'))
for RNAseqFile in RNAseqFiles:
    cancerName = os.path.basename(RNAseqFile).split('_')[0]
    cancerDict[cancerName]['RNA'] = RNAseqFile

cancerList = sorted(cancerDict.keys())

pool = Pool(processes=args.cpu)
resultList = []
miRNAlist = sorted(pairDict.keys())
for miRNAid in miRNAlist:
    geneIDlist = sorted(set(pairDict[miRNAid]))
    #result = runCoExp(miRNAid, geneIDlist)
    result = pool.apply_async(runCoExp, args=(miRNAid, geneIDlist,))
    resultList.append(result)

pool.close()
pool.join()

with open(args.output, 'w') as out:
    headerList = ['miRNAid', 'miRNAname', 'coExpGeneID', 'coExpGeneName']
    headerList.extend(list(map(lambda x:'Pan' + x, cancerList)))
    out.write('\t'.join(headerList) + '\n')
    for result in resultList:
        miRNAid, coExpGeneDict, coExpDict = result.get()
        miRNAname = miRNAnameDict[miRNAid]
        for coExpGeneID in sorted(coExpGeneDict.keys()):
            coExpGeneName = coExpGeneDict[coExpGeneID]
            coExpResultList = [miRNAid, miRNAname, coExpGeneID, coExpGeneName]
            for cancer in cancerList:
                coExpResultList.append(','.join(coExpDict[coExpGeneID][cancer]))
            out.write('\t'.join(coExpResultList) + '\n')

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))

