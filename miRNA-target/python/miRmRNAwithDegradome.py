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
                    help='The degradome-cleavage bed')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='The integrated output file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def intersect(siteDegradomeDcit, siteID, degradomeID):
    if siteID not in siteDegradomeDcit:
        siteDegradomeDcit[siteID] = [degradomeID]
    else:
        tempList = siteDegradomeDcit[siteID]
        tempList.append(degradomeID)
        siteDegradomeDcit[siteID] = tempList

starttime = datetime.datetime.now()
sys.stderr.write("Program starting!\n")

sys.stderr.write("Intersecting with degradome-seq data!\n")
pool = Pool(processes=args.cpu)
siteDegradomeDcit = Manager().dict()
command = 'bedtools intersect -a "{0}" -b "{1}" -s -wa -wb'
command = command.format(args.input, args.bed)
interResult = bytes.decode(subprocess.check_output(command, shell=True))
interList = interResult.strip().split('\n')
if len(interList) > 1:
    for inter in interList:
        row = inter.split('\t')
        siteID = row[3]
        degradomeID = row[-3]
        pool.apply_async(intersect, args=(siteDegradomeDcit, siteID, degradomeID))
    pool.close()
    pool.join()

sys.stderr.write("Writing data to ouput!\n")
with open(args.output, 'w') as out:
    with open(args.input, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            siteID = row[3]
            if siteID in siteDegradomeDcit:
                clipIDList = sorted(siteDegradomeDcit[siteID])
                degradomeIDcat = ','.join(clipIDList)
                clipNum = str(len(clipIDList))
                experimentNum = str(len(set(list(map(lambda x:x.split('-')[0], clipIDList)))))
            else:
                degradomeIDcat = 'NA'
                degradomeSiteNum = '0'
                experimentNum = '0'
            row.append(experimentNum)
            row.append(degradomeSiteNum)
            row.append(degradomeIDcat)
            out.write('\t'.join(row) + '\n')

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
