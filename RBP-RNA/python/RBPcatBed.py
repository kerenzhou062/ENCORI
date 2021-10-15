#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
import datetime
import subprocess
import multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-input', action='store', type=str,
                    help='The input bed folder')
parser.add_argument('-agoRef', action='store', type=str,
                    help='The AGO reference file')
parser.add_argument('-meta', action='store', type=str,
                    help='metadata matrix of CLIP-seq')
parser.add_argument('-RBP', nargs='*', type=str,
                    help='Designated RBP')
parser.add_argument('-test', action='store_true',
                    help='run program with test mode')
parser.add_argument('-output', action='store', type=str,
                    default="./", help='The output merge folder')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def MergeSort(bedFileList, output):
    tempList = list(map(lambda x:'"'+x+'"', bedFileList))
    bedCat = ' '.join(tempList)
    # sorted by: chr,strand,start,end
    if args.test:
        catSortCommand = 'touch {1}'.format(bedCat, output)
    else:
        catSortCommand = 'cat {0} | sort -t $\'\t\' -k1,1 -k6,6 -k2,2n -k3,3n > {1}'.format(bedCat, output)
    subprocess.run(catSortCommand, shell=True)


def MergeBed(RBP, bedFileList, bedTypeList):
    bedTypeFileList = list()
    bedgraphTypeFileList = list()
    for i in range(len(bedFileList)):
        if bedTypeList[i] == 'bed':
            bedTypeFileList.append(bedFileList[i])
        else:
            bedgraphTypeFileList.append(bedFileList[i])
    if bedTypeFileList:
        output = os.path.join(args.output, RBP + '.bed')
        MergeSort(bedTypeFileList, output)
    if bedgraphTypeFileList:
        output = os.path.join(args.output, RBP + '.bedgraph')
        MergeSort(bedgraphTypeFileList, output)


starttime = datetime.datetime.now()
sys.stderr.write("Program starting!\n")

agoRefDict = defaultdict(int)
if os.path.getsize(args.agoRef):
    with open(args.agoRef, 'r') as f:
        for line in f:
            row = line.strip()
            agoRefDict[row] = 1

findComand = 'find "{0}" -type f -name "*.bed"'.format(args.input)
findResult = bytes.decode(subprocess.check_output(findComand, shell=True))

bedList = list(map(os.path.realpath, findResult.rstrip().split('\n')))
bedDict = defaultdict(dict)
for bed in bedList:
    basename = os.path.basename(bed)
    bedDict[basename] = bed

metaDict = defaultdict(list)

with open(args.meta, 'r') as f:
    # dataSetId,Species,SeqType,GeneSymbol,Cell/Tissue,Treatment,Source,
    # Accession-GSE,Accession,FileName-RBP,Assembly,newFileName,
    # CitationForShort,Citation,PubMedID,Title
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        datasetID = row[0]
        RBP = row[3]
        fileName = row[9]
        if RBP not in agoRefDict:
            metaDict[RBP].append(bedDict[fileName])

pool = multiprocessing.Pool(processes=args.cpu)

if args.RBP:
    RBPlist = sorted(args.RBP)
else:
    RBPlist = sorted(metaDict.keys())

for RBP in RBPlist:
    bedFileList = metaDict[RBP]
    bedTypeList = list()
    for bedFile in bedFileList:
       with open(bedFile, 'r') as f:
            row = f.readline().strip().split('\t')
            if row[5] != '.':
                bedTypeList.append('bed')
            else:
                bedTypeList.append('bedgraph')
    if args.test:
        MergeBed(RBP, bedFileList, bedTypeList)
    else:
        pool.apply_async(MergeBed, (RBP, bedFileList, bedTypeList))

pool.close()
pool.join()

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
