#!/usr/bin/env python3
import os
import sys
import argparse
import re
import math
from collections import defaultdict
import datetime

parser = argparse.ArgumentParser()
parser.add_argument('-clip', action='store', type=str,
                    help='The mirTargets-clip integrate file')
parser.add_argument('-degradome', action='store', type=str,
                    help='The mirTargets-degradome integrate file')
parser.add_argument('-output', action='store', type=str, required=True,
                    help='The integrated output file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def ReconAlignment (mirnaSeq, alignment, targetSeq, watsonCrickDict, wobblePairDict):
    ## to make alignment completed, For example
    ##ugccugccgguucuaUACCUCa      *  *ugccugccgguucuaUACCUCa
    ##               ||||||       *->* :|    |::::   |||||| 
    ##uugauauguuggaugAUGGAGu      *  *uugauauguuggaugAUGGAGu
    miseqList = list(mirnaSeq)
    tarseqList = list(targetSeq)
    pairList = list(alignment)
    ## excluding the first nucleotide of the miRNA, which is not available for pairing
    ## ref. Metazoan MicroRNAs. Cell, 2018
    for i in range(len(miseqList) - 1):
        mread = miseqList[i]
        tread = tarseqList[i]
        pair = pairList[i]
        if pair != '|' and tread in watsonCrickDict[mread]:
            pairList[i] = '|'
        elif pair != '|' and tread in wobblePairDict[mread]:
            pairList[i] = ':'
    ## set the first nucloetide to be unpaired
    pairList[-1] = ' '
    reconAlign = ''.join(pairList)
    return reconAlign

def CreatWatsonCrickDict (watsonList, crickList):
    if len(watsonList) != len(crickList):
        sys.stderr.write("watson and crick list should be in the same length")
        sys.exit()
    pairDict = defaultdict(set)
    for i in range(len(watsonList)):
        wlist = [watsonList[i].upper(), watsonList[i].lower()]
        clist = [crickList[i].upper(), crickList[i].lower()]
        for j in range(len(wlist)):
            for k in range(len(clist)):
                pairDict[wlist[j]].update([clist[k]])
                pairDict[clist[k]].update([wlist[j]])
    return pairDict

def MirnaPosTag(mirnaSeq, pos):
    mirnaSeqRev = mirnaSeq[::-1]
    indexFlag = 0
    count = 0
    for i in range(len(mirnaSeqRev)):
        if count == pos:
            indexFlag = i
            break
        else:
            if mirnaSeqRev[i] == '-':
                continue
            else:
                count += 1
    return indexFlag


def MirTargetSitesType(mirnaSeq, alignment, targetSeq):
    ## ref. Metazoan MicroRNAs, Cell, 2018
    ###         3' UTR
    ### NNNNNN.. Offset 6mer
    ### .NNNNNN. 6mer
    ### .NNNNNNA 7mer-A1
    ### NNNNNNN. 7mer-m8
    ### NNNNNNNA 8mer
    ### NNNNNNNN (miRNA)
    mirnaSeqRev = mirnaSeq[::-1]
    alignRev = alignment[::-1]
    targetSeqRev = targetSeq[::-1].upper()
    ## seed region, 2-7 nucleotide of miRNA
    seedMatchCount = alignRev[1:7].count('|')
    ## the last nucloetide of target sequence
    lastAFlag = False
    if targetSeqRev[0] == 'A':
        lastAFlag = True
    ## whether the 8th nucleotide of miRNA matched to target
    match8thFlag = False
    if alignRev[7] == '|':
        match8thFlag = True
    offsetSeedMatchCount = alignRev[2:7].count('|')
    ## determin the target sites
    if seedMatchCount == 6:
        ## canonical sites
        if lastAFlag is False and match8thFlag is False:
            siteType = '6mer'
        elif lastAFlag is True and match8thFlag is False:
            siteType = '7mer-A1'
        elif lastAFlag is False and match8thFlag is True:
            siteType = '7mer-m8'
        elif lastAFlag is True and match8thFlag is True:
            siteType = '8mer'
        else:
            siteType = 'other'
    elif seedMatchCount == 5:
        ## for noncanonical sites
        ## only 1 mismatch in seed region is allowed
        ##1 Productive 3′-supplementary pairing typically centers on nucleotides 13–16
        index13Flag = MirnaPosTag(mirnaSeq, 13)
        prime3CenterMatch = alignRev[index13Flag:index13Flag + 4].count('|')
        ##2 A centered site is one that lacks perfect seed pairing and 3'-compensatory pairing but instead has 11-12 contiguous Watson-Crick pairs to miRNA positions 4-15
        index4Flag = MirnaPosTag(mirnaSeq, 4)
        index15Flag = MirnaPosTag(mirnaSeq, 15)
        centerMatchCount = alignRev[index4Flag:index15Flag + 1].count('|')
        if prime3CenterMatch == 4 or (centerMatchCount == 11 or centerMatchCount == 12):
            siteType = 'noncanonical'
        else:
            siteType = 'non-seed'
    else:
        ## canonical sites, offset-6mer
        if offsetSeedMatchCount == 5 and match8thFlag is True:
            siteType = 'offset-6mer'
        else:
            siteType = 'non-seed'
    return siteType

def FindBulge(alignment):
    #alignment:3'->5', seed region, 2nd~7th
    alignRev = alignment[::-1].replace(':', '|')
    matchIter = re.finditer(r'\|+', alignRev)
    matchIndexList = [[m.start(), m.end() - 1] for m in matchIter]
    seedEnd = 6
    bulgeStart = 0
    bulgeStartPos = 0
    for i in range(len(matchIndexList)):
        mstart = matchIndexList[i][0]
        mend = matchIndexList[i][1]
        if mstart <= seedEnd and seedEnd <= mend:
            bulgeStart = mend + 1
            bulgeStartPos = i
    if len(matchIndexList) > (bulgeStartPos + 1):
        bulgeEnd = matchIndexList[bulgeStartPos + 1][0] - 1
    else:
        bulgeEnd = len(alignment) - 1
    ## bulgeStart can be equal to bulgeEnd， like [7, 7]
    bulgeRegion = [bulgeStart, bulgeEnd]
    return bulgeRegion

def TagRegionByBulge(alignRev, regex, bulgeRegion):
    regexIter = re.finditer(r'{0}'.format(regex), alignRev)
    regexIndexList = [[m.start(), m.end() - 1] for m in regexIter]
    bulgeStart = bulgeRegion[0]
    bulgeEnd = bulgeRegion[1]
    bulgeSeedEnd = bulgeStart + 4
    extraIndexList = []
    ## add tag to the region
    ## "front" means the region located in front of the bulge region
    ## "in" means the region located in the seed bulge region (<= 6nt)
    ## "beyond" means the region located in bulge region but beyond the seed bulge region
    ## "behind" means the region located behind the bulge region
    for i in range(len(regexIndexList)):
        regionStart = regexIndexList[i][0]
        regionEnd = regexIndexList[i][1]
        if regionEnd < bulgeStart:
            regexIndexList[i].append('front')
        elif regionStart > bulgeEnd:
            regexIndexList[i].append('behind')
        else:
            ## this should be mismatch or gap, no match in the bulge region
            if regionEnd <= bulgeSeedEnd:
                regexIndexList[i].append('in')
            elif regionStart > bulgeSeedEnd:
                regexIndexList[i].append('beyond')
            else:
                ##regionStart <= bulgeSeedEnd and regionEnd > bulgeSeedEnd
                regexIndexList[i] = [regionStart, bulgeSeedEnd, 'in']
                extraIndexList.append([bulgeSeedEnd + 1, regionEnd, 'beyond'])
    regexIndexList.extend(extraIndexList)
    return regexIndexList

def RegionClassifier(mirnaSeq, alignment, targetSeq):
    mirnaSeqList = list(mirnaSeq[::-1])
    pairList = list(alignment[::-1])
    targetSeqList = list(targetSeq[::-1])
    ## to identify mismatch and gap from unpairs
    for i in range(len(pairList)):
        mirnaNuc = mirnaSeqList[i]
        pair = pairList[i]
        targetNuc = targetSeqList[i]
        if pair == ' ':
            if mirnaNuc == '-' or targetNuc == '-':
                pairList[i] = '-'
    alignRev = ''.join(pairList)
    ## classify and label the region
    bulgeRegion = FindBulge(alignment)
    matchTagList = TagRegionByBulge(alignRev, '\|+', bulgeRegion)
    wobbleTagList = TagRegionByBulge(alignRev, ':+', bulgeRegion)
    mismatchTagList = TagRegionByBulge(alignRev, '\s+', bulgeRegion)
    gapTagList = TagRegionByBulge(alignRev, '-+', bulgeRegion)
    ## store the labeled region
    regionDict = {}
    regionDict['match'] = matchTagList
    regionDict['wobble'] =wobbleTagList
    regionDict['mismatch'] = mismatchTagList
    regionDict['gap'] = gapTagList
    return regionDict

def CalculateTdmdScore(mirnaSeq, alignment, targetSeq):
    # calculate a score for TDMD
    ## mirnaSeq: 3'->5', alignment:3'->5', targetSeq:5'->3'
    ## like:UUGUUGUUUUAGUG-AU-CAGAAGGU,||||||||||||||    ||||||| ,AACAACAAAAUCACCAAUGUCUUCCA
    ## score at 3' region: eg. match_score = continousMatchCount * matchScore
    ## construct penalty dictionary
    penaltyDict = {}
    penaltyDict['match'] = 2
    penaltyDict['wobble'] = 1
    penaltyDict['gap'] = -2
    penaltyDict['mismatch'] = -1
    reAlignment = ReconAlignment(mirnaSeq, alignment, targetSeq, watsonCrickDict, wobblePairDict)
    siteType = MirTargetSitesType(mirnaSeq, reAlignment, targetSeq)
    if (siteType == "non-seed"):
        startScore = -2 * len(targetSeq)
    else:
        startScore = 0
    regionDict = RegionClassifier(mirnaSeq, alignment, targetSeq)
    ## {'match': [[1, 7, 'front'], [12, 25, 'behind']], 'wobble': [], 'mismatch': [[0, 0, 'front'], [9, 10, 'in']], 'gap': [[8, 8, 'in'], [11, 11, 'in']]}
    ## calculate the penalty score
    score = 0
    for regionType in sorted(regionDict.keys()):
        if bool(regionDict[regionType]):
            for region in regionDict[regionType]:
                ## ignore the fisrt nucleotide of miRNA
                if regionType == 'mismatch':
                    if region[0] == 0 and region[1] == 0:
                        continue
                    elif region[0] == 0 and region[1] > 0:
                        region[0] = 1
                start = region[0]
                end = region[1]
                regionTag = region[2]
                contigousCount = 1
                for i in range(end - start + 1):
                    if regionTag == 'beyond':
                        score += 2 * penaltyDict[regionType]
                    if regionTag == 'behind':
                        score += math.sqrt(contigousCount) * penaltyDict[regionType]
                        contigousCount += 1
                    else:
                        score += penaltyDict[regionType]
    ## eliminate the 1st nucleotide
    score = score - penaltyDict['mismatch'] + startScore
    ## calculate the length of miRNA
    mirnaLength = 0
    for i in range(len(mirnaSeq)):
        if mirnaSeq[i] != '-':
            mirnaLength += 1
    ## normalize the score by miRNA length
    tdmdScore = round(score / mirnaLength, 4)
    return tdmdScore

'''
example for tdmdScore
eg1: miR-7 and Cyrano
mirnaSeq =  'UUGUUGUUUUAGUG-AU-CAGAAGGU'
alignment = '||||||||||||||    ||||||| '
targetSeq = 'AACAACAAAAUCACCAAUGUCUUCCA'
CalculateTdmdScore(mirnaSeq, alignment, targetSeq)
tdmd score: 3.4247

eg2: miR-29b and NREP
mirnaSeq =  'GUAACCGAU---GGAUGGUGCUA'
alignment = ':  || |||   ::|||||||||'
targetSeq = 'UUGUGACUAAAGUUUACcacgaU'
CalculateTdmdScore(mirnaSeq, alignment, targetSeq)
tdmd score: 1.2853

eg3: ramdon pairing
mirnaSeq =  'CCAGUGUUGGCAUG-AU-CAGAAGGU'
alignment = '    ||||    ||    ||||||| '
targetSeq = 'AACAACAAAAUCACCAAUGUCUUCCA'
CalculateTdmdScore(mirnaSeq, alignment, targetSeq)
tdmd score: 0.5762
'''

starttime = datetime.datetime.now()
sys.stderr.write("Program starting!\n")

sys.stderr.write("reading mirTargets-degradome integrate file!\n")

with open(args.degradome, 'r') as f:
    degradList = f.readlines()

sys.stderr.write("reading mirTargets-clip integrate file!\n")
with open(args.clip, 'r') as f:
    clipList = f.readlines()

sys.stderr.write("Writing data to ouput!\n")

watsonList = ['A', 'U', 'C', '-']
crickList = ['T', 'A', 'G', '-']
watsonCrickDict = CreatWatsonCrickDict(watsonList, crickList)

wb1List = ['G']
wb2List = ['U']
wobblePairDict = CreatWatsonCrickDict(wb1List, wb2List)

with open(args.output, 'w') as out:
    for i in range(len(clipList)):
        clipRow = clipList[i].strip().split('\t')
        degradomeRow = degradList[i].strip().split('\t')
        row = clipRow[0:11]
        mirnaSeq = row[8]
        alignment = row[9]
        targetSeq = row[10]
        tdmdScore = CalculateTdmdScore(mirnaSeq, alignment, targetSeq)
        clipExpNum = clipRow[11]
        clipSiteNum = clipRow[12]
        clipIDcat = clipRow[13]
        degradExpNum = degradomeRow[11]
        degradSiteNum = degradomeRow[12]
        degradIDcat = degradomeRow[13]
        row.extend([str(tdmdScore), clipExpNum, clipSiteNum, clipIDcat, degradExpNum, degradSiteNum, degradIDcat])
        out.write('\t'.join(row) + '\n')

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!\n")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
