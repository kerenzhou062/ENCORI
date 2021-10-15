#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
import subprocess
from multiprocessing import Pool, Manager
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-input', action='store', type=str,
                     help='miRNA binding-target file')
parser.add_argument('-bed', action='store', type=str,
                    help='The folder contain all processed-CLIP-bed files')
parser.add_argument('-meta', action='store', type=str,
                    help='metadata matrix of CLIP-seq')
parser.add_argument('-agoRef', action='store', type=str, required=True,
                    help='The reference file for AGO protein')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='The integrated output file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def intersect(siteClipDcit, bedA, bedB, strand):
    command = ''
    if strand != '.':
        command = 'bedtools intersect -a "{0}" -b "{1}" -s -wa -wb'
    else:
        command = 'bedtools intersect -a "{0}" -b "{1}" -wa -wb'
    command = command.format(bedA, bedB)
    interResult = bytes.decode(subprocess.check_output(command, shell=True))
    interList = interResult.strip().split('\n')
    if len(interList) > 1:
        for inter in interList:
            row = inter.split('\t')
            siteID = row[3]
            clipID = row[-3]
            if siteID not in siteClipDcit:
                siteClipDcit[siteID] = [clipID]
            else:
                tempList = siteClipDcit[siteID]
                tempList.append(clipID)
                siteClipDcit[siteID] = tempList

starttime = datetime.datetime.now()
sys.stderr.write("Program starting!\n")

agoDict = defaultdict(int)
with open(args.agoRef, 'r') as f:
    for line in f:
        agoDict[line.strip()] = 1

findComand = 'find "{0}" -type f -name "*.bed"'.format(args.bed)
findResult = bytes.decode(subprocess.check_output(findComand, shell=True))

bedList = list(map(os.path.realpath, findResult.rstrip().split('\n')))
bedDict = defaultdict(dict)
for bed in bedList:
    basename = os.path.basename(bed)
    bedDict[basename] = bed

sys.stderr.write("Intersecting with AGO-CLIP-seq data!\n")

pool = Pool(processes=args.cpu)
siteClipDcit = Manager().dict()
with open(args.meta, 'r') as f:
    #dataSetId,Species,SeqType,GeneSymbol,Cell/Tissue,Treatment,Source,Accession-GSE,Accession,
    #FileName-RBP,Assembly,newFileName,CitationForShort,Citation,PubMedID,Title
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        RBP = row[3]
        if RBP not in agoDict:
            continue
        bedName = row[9]
        clipBed = bedDict[bedName]
        with open(clipBed, 'r') as f:
            strand = f.readline().strip().split('\t')[-1]
        #intersect(args.input, clipBed, strand)
        pool.apply_async(intersect, args=(siteClipDcit, args.input, clipBed, strand))
pool.close()
pool.join()

sys.stderr.write("Writing data to ouput!\n")
with open(args.output, 'w') as out:
    with open(args.input, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            siteID = row[3]
            if siteID in siteClipDcit:
                clipIDList = sorted(siteClipDcit[siteID])
                clipIDcat = ','.join(clipIDList)
                clipSiteNum = str(len(clipIDList))
                experimentNum = str(len(set(list(map(lambda x:x.split('-')[0], clipIDList)))))
            else:
                clipIDcat = 'NA'
                clipSiteNum = '0'
                experimentNum = '0'
            row.append(experimentNum)
            row.append(clipSiteNum)
            row.append(clipIDcat)
            out.write('\t'.join(row) + '\n')

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
