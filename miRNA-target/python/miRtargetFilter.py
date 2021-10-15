#!/usr/bin/env python3
import os
import sys
import argparse
import tempfile
from collections import defaultdict
import subprocess
from multiBioPro import BedMan

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The input file')
parser.add_argument('-fasta', action='store', type=str,
                    help='The fasta of genome')
parser.add_argument('-output', action='store', type=str,
                    default="output.txt", help='The output file')

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

watsonList = ['A', 'U', 'C', '-']
crickList = ['T', 'A', 'G', '-']
watsonCrickDict = CreatWatsonCrickDict(watsonList, crickList)

wb1List = ['G']
wb2List = ['U']
wobblePairDict = CreatWatsonCrickDict(wb1List, wb2List)

# get target site sequence
getFastaCommand = "bedtools getfasta -bed {0} -fi {1} -s -bedOut 2> /dev/null".format(args.input, args.fasta);
resultLineList = bytes.decode(subprocess.check_output(getFastaCommand, shell=True)).split('\n')
miRdict = defaultdict(str)
with open (args.output, 'w') as out:
    for line in resultLineList:
        if line == '':
            continue
        row = line.strip().split('\t')
        ## test whether the sequence of target bed matched to genome sequence
        targetSeq = row[10].replace('-', '').lower().replace('t', 'u')
        targetFasta = row[-1].lower().replace('t', 'u')
        if targetFasta in targetSeq:
            pass
        elif targetSeq in targetFasta:
            pass
        else:
            continue
        ## get the target site type
        ## seed region, 2:7
        mirnaSeq = row[8]
        alignment = row[9]
        targetSeq = row[10]
        reAlignment = ReconAlignment(mirnaSeq, alignment, targetSeq, watsonCrickDict, wobblePairDict)
        row[4] = MirTargetSitesType(mirnaSeq, reAlignment, targetSeq)
        ## discard non-seed miRNA-target
        if row[4] == 'non-seed':
            continue
        tempRow = row[0:3]
        tempRow.extend(row[4:11])
        ID = row[3]
        tempKey = '\t'.join(tempRow)
        if tempKey not in miRdict:
            miRdict[tempKey] = ID
            row[9] = reAlignment
            row[8] = row[8].replace('T', 'U')
            row[8] = row[8].replace('t', 'u')
            row[10] = row[10].replace('T', 'U')
            row[10] = row[10].replace('t', 'u')
            out.write('\t'.join(row[0:11]) + '\n')
