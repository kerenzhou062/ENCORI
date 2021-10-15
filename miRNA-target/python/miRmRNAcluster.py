#!/usr/bin/env python3
import os
import sys
import argparse
from collections import defaultdict
from functools import reduce
import datetime
from multiBioPro import BedMan

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The integrated mirTarget-seq file')
parser.add_argument('-prefix', action='store', type=str,
                    help='The prefix for cluster ID')
parser.add_argument('-meta', action='store', type=str,
                    help='metadata matrix of CLIP-seq')
parser.add_argument('-ref', action='store', type=str,
                    help='The reference file contain binding ids and predict \
                    software')
parser.add_argument('-output', action='store', type=str,
                    default="cluster.txt", help='The output cluster file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def DictWrite(miRNA, cluster, tdmdScoreList, siteList, siteIDlist, siteClipDict, siteDegraDict, Object):
    # determin the maximum tdmdScore
    tdmdScore = max(tdmdScoreList)
    #transform siteClipDict to RBP and corresponding clipID(peak)
    #statistics for software peak number
    allClipIDlist = list()
    siteIDclipList = [[] for i in range(softwareNum)]
    RBPdict = defaultdict(list)
    for siteID in siteClipDict:
        clipIDcat = siteClipDict[siteID]
        if clipIDcat == 'NA':
            continue
        clipIDlist = clipIDcat.split(',')
        software = refDict[siteID]
        index = softwares.index(software)
        siteIDclipList[index].extend(clipIDlist)
        allClipIDlist.extend(clipIDlist)
        for clipID in clipIDlist:
            datasetID = clipID.split('-')[0]
            RBP = metaDict[datasetID]
            RBPdict[RBP].append(clipID)
    allClipIDnum = len(set(allClipIDlist))
    allClipExpNum = len(set(map(lambda x:x.split('-')[0], allClipIDlist)))
    RBPList = sorted(RBPdict.keys())
    if RBPList:
        RBPclipList = list()
        RBPClipNumList = list()
        for RBP in RBPList:
            clipIDlist = sorted(set(RBPdict[RBP]))
            clipIDlistCat = ':'.join(clipIDlist)
            RBPclipList.append(clipIDlistCat)
            RBPClipNumList.append(str(len(clipIDlist)))
        RBPcat = ','.join(RBPList)
        RBPClipNumCat = ','.join(RBPClipNumList)
        RBPclipIDcat = ','.join(RBPclipList)
    else:
        RBPcat = 'NA'
        RBPClipNumCat = 'NA'
        RBPclipIDcat = 'NA'
    allDegraIDlist = list()
    siteDegraList = [[] for i in range(softwareNum)]
    for siteID in siteDegraDict:
        degraIDcat = siteDegraDict[siteID]
        if degraIDcat == 'NA':
            continue
        degraIDlist = degraIDcat.split(',')
        allDegraIDlist.extend(degraIDlist)
        software = refDict[siteID]
        index = softwares.index(software)
        siteDegraList[index].extend(degraIDlist)
    if allDegraIDlist:
        uniqDegraIDList = sorted(set(allDegraIDlist))
        degraIDcat = ','.join(uniqDegraIDList)
        degraExpNum = len(set(map(lambda x:x.split('-')[0], uniqDegraIDList)))
        degraIDNum = len(uniqDegraIDList)
    else:
        degraIDcat = 'NA'
        degraExpNum = 0
        degraIDNum = 0

    siteWithSeqList = ['' for i in range(softwareNum)]
    for i in range(softwareNum):
        if siteIDclipList[i]:
            clipExpnum = str(len(set(map(lambda x:x.split('-')[0], siteIDclipList[i]))))
        else:
            clipExpnum = '0'
        if siteDegraList[i]:
            degraExpnum = str(len(set(map(lambda x:x.split('-')[0], siteDegraList[i]))))
        else:
            degraExpnum = '0'
        siteWithSeqList[i] = ','.join([clipExpnum, degraExpnum])

    clusterID = args.prefix + '%09d' % cluster
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
    #
    siteIDtempList = [[] for i in range(softwareNum)]
    for siteID in siteIDlist:
        software = refDict[siteID]
        index = softwares.index(software)
        siteIDtempList[index].append(siteID)
    siteIDcoutList = list()
    siteIDcatList = list()
    softwareCount = 0
    for siteIDtemp in siteIDtempList:
        if siteIDtemp:
            softwareCount += 1
            siteIDcat = ','.join(siteIDtemp)
            siteIDcount = len(siteIDtemp)
        else:
            siteIDcat = 'NA'
            siteIDcount = 0
        siteIDcoutList.append(siteIDcount)
        siteIDcatList.append(siteIDcat)
    row = [chro]
    row.extend(intersectCluster)
    row.extend([clusterID, softwareCount, strand, miRNA ,name])
    row.extend(mergeCluster)
    row.append(tdmdScore)
    row.extend(siteIDcoutList)
    row.extend(siteIDcatList)
    row.extend(siteWithSeqList)
    row.extend([allClipExpNum, allClipIDnum, RBPcat, RBPClipNumCat,
        RBPclipIDcat, degraExpNum, degraIDNum, degraIDcat])
    Object.write('\t'.join(list(map(str,row))) + '\n')

starttime = datetime.datetime.now()
sys.stderr.write("program starting!\n")

refDict = defaultdict(dict)
softwares = list()
with open(args.ref, 'r') as f:
    for line in f.readlines():
        row = line.strip().split('\t')
        refDict[row[0]] = row[1]
        if row[1] not in softwares:
            softwares.append(row[1])

metaDict = defaultdict(str)
with open(args.meta, 'r') as f:
    # dataSetId,Species,SeqType,GeneSymbol,Cell/Tissue,Treatment,Source,
    # Accession-GSE,Accession,FileName-RBP,Assembly,newFileName,
    # CitationForShort,Citation,PubMedID,Title
    __ = f.readline()
    for line in f:
        row = line.strip().split('\t')
        datasetID = row[0]
        RBP = row[3]
        metaDict[datasetID] = RBP

softwares = sorted(softwares)
softwareNum = len(softwares)

mirIDnameDict = defaultdict(dict)
intersectDict = defaultdict(dict)
cluster = 1
with open(args.output, 'w') as out:
    headerList = ['chromosome', 'narrowStart', 'narrowEnd', 'clusterID',
        'softwareCount', 'strand', 'miRNAid', 'miRNAname', 'broadStart',
        'broadEnd','tdmdScore']
    headerList.extend(softwares)
    headerList.extend(list(map(lambda x:x+'ID', softwares)))
    headerList.extend(list(map(lambda x:x+'SeqStats', softwares)))
    headerList.extend(['clipExpNum', 'clipIDnum', 'RBP', 'RBPclipNum',
        'RBPclipID', 'degraExpNum', 'degraIDNum', 'degraID'])
    out.write('\t'.join(headerList) + '\n')
    with open(args.input, 'r') as f:
        #chr1:879743:879765:HURNATH00002260 0 :-:MIMAT0000062:hsa-let-7a-5p:
        #uugauauguuGGAUGAUGGAGu::| || |||||:-------ccuUCCACCACCUCu:3.383:
        #4:6:SBDH206-4695,...:5:36:DEHSA0006-33068,...
        last = f.readline().strip().split('\t')
        last[1] = int(last[1])
        last[2] = int(last[2])
        tdmdScore = float(last[11])
        tdmdScoreList = []
        siteList = []
        siteIDlist = []
        tdmdScoreList.append(tdmdScore)
        siteList.append([last[0], last[1], last[2], last[5]])
        siteIDlist.append(last[3])
        mirIDnameDict[last[6]] = last[7]
        siteClipDict = {last[3]:last[14]}
        siteDegraDict = {last[3]:last[17]}
        intersectDict[cluster] = False

        for line in f:
            row = line.strip().split('\t')
            row[1] = int(row[1])
            row[2] = int(row[2])
            chro = row[0]
            start = row[1]
            end = row[2]
            siteID = row[3]
            strand = row[5]
            miRNA = row[6]
            name = row[7]
            tdmdScore = float(row[11])

            if miRNA != last[6]:
                DictWrite(miRNA, cluster, tdmdScoreList, siteList, siteIDlist, siteClipDict, siteDegraDict, out)
                cluster += 1
                tdmdScoreList = []
                siteList = []
                siteIDlist = []
                tdmdScoreList.append(tdmdScore)
                siteList.append([chro, start, end, strand])
                siteIDlist.append(siteID)
                siteClipDict = {siteID:row[14]}
                siteDegraDict = {siteID:row[17]}
                intersectDict[cluster] = False
                last = row
            else:
                #if same chromosome and strand
                if last[0] == chro and last[5] == strand:
                    intersect = BedMan.intersect([start, end], [last[1], last[2]])
                    # if intersect with last line
                    if intersect:
                        # if common-cluster-interval exits
                        if intersectDict[cluster]:
                            tempIntersect = BedMan.intersect(intersect, intersectDict[cluster])
                            # if the instersect overlap with common-cluster-interval
                            if tempIntersect:
                                tdmdScoreList.append(tdmdScore)
                                siteList.append([chro, start, end, strand])
                                siteIDlist.append(siteID)
                                siteClipDict[siteID] = row[14]
                                siteDegraDict[siteID] = row[17]
                                intersectDict[cluster] = tempIntersect
                                last = row
                            else:
                                DictWrite(miRNA, cluster, tdmdScoreList, siteList, siteIDlist, siteClipDict, siteDegraDict, out)
                                cluster += 1
                                tdmdScoreList = []
                                siteList = []
                                siteIDlist = []
                                tdmdScoreList.append(tdmdScore)
                                siteList.append([chro, start, end, strand])
                                siteIDlist.append(siteID)
                                siteClipDict = {siteID:row[14]}
                                siteDegraDict = {siteID:row[17]}
                                intersectDict[cluster] = False
                                last = row
                        else:
                            tdmdScoreList.append(tdmdScore)
                            siteList.append([chro, start, end, strand])
                            siteIDlist.append(siteID)
                            siteClipDict[siteID] = row[14]
                            siteDegraDict[siteID] = row[17]
                            intersectDict[cluster] = intersect
                            last = row
                    else:
                        DictWrite(miRNA, cluster, tdmdScoreList, siteList, siteIDlist, siteClipDict, siteDegraDict, out)
                        cluster += 1
                        tdmdScoreList = []
                        siteList = []
                        siteIDlist = []
                        tdmdScoreList.append(tdmdScore)
                        siteList.append([chro, start, end, strand])
                        siteIDlist.append(siteID)
                        siteClipDict = {siteID:row[14]}
                        siteDegraDict = {siteID:row[17]}
                        intersectDict[cluster] = False
                        last = row
                else:
                    DictWrite(miRNA, cluster, tdmdScoreList, siteList, siteIDlist, siteClipDict, siteDegraDict, out)
                    cluster += 1
                    tdmdScoreList = []
                    siteList = []
                    siteIDlist = []
                    tdmdScoreList.append(tdmdScore)
                    siteList.append([chro, start, end, strand])
                    siteIDlist.append(siteID)
                    siteClipDict = {siteID:row[14]}
                    siteDegraDict = {siteID:row[17]}
                    intersectDict[cluster] = False
                    last = row

            if miRNA not in mirIDnameDict:
                mirIDnameDict[miRNA] = name

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
