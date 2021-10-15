#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
import datetime
import subprocess
from glob import glob
from functools import reduce
from multiprocessing import Pool, Manager
from multiBioPro import BedMan

parser = argparse.ArgumentParser()
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-distance', action='store', type=int,
                    default=20, help='-d for bedtools merge')
parser.add_argument('-input', action='store', type=str,
                    help='The input bed-merge folder')
parser.add_argument('-prefix', action='store', type=str,
                    help='The prefix for cluster ID')
parser.add_argument('-test', action='store_true',
                    help='run program with test mode')
parser.add_argument('-output', action='store', type=str,
                    default="./", help='The output merge folder')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def DictWrite(RBP, prefix, output, cluster, siteList, siteIDlist):
    locusList = list()
    for site in siteList:
        locusList.append([site[1], site[2]])
    if len(locusList) > 1:
        mergeCluster = reduce(BedMan.merge, locusList)
        intersectCluster = reduce(BedMan.intersect, locusList)
    else:
        mergeCluster = locusList[0]
        intersectCluster = locusList[0]
    chro = siteList[0][0]
    strand = siteList[0][-1]
    siteIDnum = len(siteIDlist)
    expNum = len(set(list(map(lambda x:x.split('-')[0], siteIDlist))))
    siteIDcat = ','.join(sorted(siteIDlist))
    clusterID = prefix + '%09d' % cluster

    outputList = [chro]
    outputList.extend(intersectCluster)
    outputList.extend([clusterID, expNum, strand])
    outputList.extend(mergeCluster)
    outputList.extend([RBP, siteIDnum, siteIDcat])
    output.write('\t'.join(list(map(str, outputList))) + '\n')

def BedCluster(RBP, bedFile, outputFile):
    intersectDict = defaultdict(dict)
    prefix = RBP + ':' + args.prefix
    ## merge peak with -d 20
    command = 'sort -t$\'\\t\' -k1,1 -k2,2n {} | bedtools merge -i stdin -d {} -s -c 4,5,6 -o collapse,sum,distinct'.format(bedFile, args.distance)
    runResList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
    mergePeakRowList = list(map(lambda x:x.split('\t'), list(filter(lambda x: bool(x), runResList))))
    with open(outputFile, 'w') as output:
        for i in range(len(mergePeakRowList)):
            #['chr1', '630790', '630792', 'SBDH1356-T-11,SBDH1356-T-12', '396', '-']
            cluster = i + 1
            mergePeakRow = mergePeakRowList[i]
            mergeChr = mergePeakRow[0]
            mergeStart = mergePeakRow[1]
            mergeEnd = mergePeakRow[2]
            siteIDlist = mergePeakRow[4].split(',')
            siteIDnum = len(siteIDlist)
            expNum = len(set(list(map(lambda x:x.split('-')[0], siteIDlist))))
            clusterID = prefix + str(cluster)
            # outputlist
            outputList = mergePeakRow[0:3]
            outputList.extend([clusterID, expNum, mergePeakRow[3]])
            outputList.extend(mergePeakRow[1:3])
            outputList.extend([RBP, siteIDnum, mergePeakRow[4]])
            output.write('\t'.join(list(map(str, outputList))) + '\n')

starttime = datetime.datetime.now()
sys.stderr.write("Program starting!\n")

bedFileList = glob(os.path.join(args.input, '*.bed'))

pool = Pool(processes=args.cpu)
for bedFile in bedFileList:
    baseName = os.path.basename(bedFile)
    RBP = os.path.splitext(baseName)[0]
    outputFile = os.path.join(args.output, RBP + '.txt')
    if args.test:
        BedCluster(RBP, bedFile, outputFile)
    else:
        pool.apply_async(BedCluster, args=(RBP, bedFile, outputFile))

pool.close()
pool.join()

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
