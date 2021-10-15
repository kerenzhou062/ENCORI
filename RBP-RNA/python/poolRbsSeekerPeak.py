#!/usr/bin/env python3
import os
import sys
import re
import argparse
from glob import glob
from collections import defaultdict
import subprocess
from multiprocessing import Pool, Manager
import datetime
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-meta', action='store', type=str,
                    help='metadata matrix')
parser.add_argument('-max', action='store', type=int,
                    default=200, help='The maxmum length of peak')
parser.add_argument('-min', action='store', type=int,
                    default=4, help='The minimum length of peak')
parser.add_argument('-bed', action='store', type=str,
                    help='folder contain bed files')
parser.add_argument('-output', action='store', type=str,
                    default="./", help='The output directory')
args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def GetBedLocusCount(bedFile):
    count = 0
    with open(bedFile, 'r') as f:
        for line in f:
            count += 1
    return count

def RebuildSiteBed(inputFile, outputFile, datasetId, clipSiteType):
    clipSiteType = clipSiteTypeShortDict[clipSiteType]
    shellCommand = "awk 'BEGIN{{FS=\"\t\";OFS=\"\t\";}}{{if ($1 !~ \"^#\"){{print $1,$2,$3,$4,$5,$6}}}}' {0} |".format(inputFile)
    shellCommand += "sort -t $'\t' -k1,1 -k2,2n | awk 'BEGIN{{OFS=\"\t\";FS=\"\t\"}}"
    shellCommand += "{{$4=\"{0}-{1}-\"FNR;print $0}}' > {2}".format(datasetId, clipSiteType, outputFile)
    __ = subprocess.check_output(shellCommand, shell=True)

def RebuildPeakBed(inputFile, outputFile, datasetId, minPeakLen, maxPeakLen, clipSiteType):
    clipSiteType = clipSiteTypeShortDict[clipSiteType]
    shellCommand = "awk 'BEGIN{{FS=\"\t\";OFS=\"\t\";}}{{if ($1 !~ \"^#\"){{print $1,$2,$3,$4,$5,$6}}}}' {0} |".format(inputFile)
    shellCommand += "sort -t $'\t' -k1,1 -k2,2n | awk 'BEGIN{{OFS=\"\t\";FS=\"\t\"}}"
    shellCommand += "{{peakLen=$3-$2;if(peakLen >= {0} && peakLen <= {1})".format(minPeakLen, maxPeakLen)
    shellCommand += "{{$4=\"{0}-{1}-\"FNR;print $0}}}}' > {2}".format(datasetId, clipSiteType, outputFile)
    __ = subprocess.check_output(shellCommand, shell=True)

def KeepNonOverlapPeak(siteBed, peakBed, outputBed):
    shellCommand = "bedtools intersect -a {0} -b {1} -s -v > {2}".format(siteBed, peakBed, outputBed)
    __ = subprocess.check_output(shellCommand, shell=True)

def CombineSitePeakBed(siteBed, peakBed, outputBed):
    if siteBed != 'none' and peakBed != 'none':
        shellCommand = "cat {0} {1} | sort -t $'\t' -k1,1 -k2,2n > {2}".format(siteBed, peakBed, outputBed)
    elif siteBed != 'none':
        shellCommand = "cat {0} | sort -t $'\t' -k1,1 -k2,2n > {1}".format(siteBed, outputBed)
    elif peakBed != 'none':
        shellCommand = "cat {0} | sort -t $'\t' -k1,1 -k2,2n > {1}".format(peakBed, outputBed)
    __ = subprocess.check_output(shellCommand, shell=True)

def PoolBed(bedInfoDict, mergeBed, minPeakLen, maxPeakLen, outputDir):
    ## bedTypeList, P|M|T
    inforDict = bedInfoDict[mergeBed]
    datasetId = inforDict['datasetId']
    if 'peak' in inforDict:
        peakFlag = True
    else:
        peakFlag = False
    if 'site' in inforDict:
        siteFlag = True
    else:
        siteFlag = False
    if siteFlag and peakFlag:
        ## for mutation|truncations|deletion bed, rename the peak_id
        clipSiteType = inforDict['site'][0]
        siteBed = inforDict['site'][1]
        formatSiteBed = tempfile.NamedTemporaryFile(delete=True)
        RebuildSiteBed(siteBed, formatSiteBed.name, datasetId, clipSiteType)
        ## remove the peaks that overlap with mutation|truncations bed
        originalPeakBed = inforDict['peak'][1]
        kepetBed = tempfile.NamedTemporaryFile(delete=True)
        KeepNonOverlapPeak(originalPeakBed, formatSiteBed.name, kepetBed.name)
        peakBed = kepetBed.name
        clipSiteType = inforDict['peak'][0]
        formatPeakBed = tempfile.NamedTemporaryFile(delete=True)
        RebuildPeakBed(peakBed, formatPeakBed.name, datasetId, minPeakLen, maxPeakLen, clipSiteType)
    elif siteFlag is False and peakFlag is True:
        clipSiteType = inforDict['peak'][0]
        peakBed = inforDict['peak'][1]
        formatPeakBed = tempfile.NamedTemporaryFile(delete=True)
        RebuildPeakBed(peakBed, formatPeakBed.name, datasetId, minPeakLen, maxPeakLen, clipSiteType)
    elif siteFlag is True and peakFlag is False:
        clipSiteType = inforDict['site'][0]
        siteBed = inforDict['site'][1]
        formatSiteBed = tempfile.NamedTemporaryFile(delete=True)
        RebuildSiteBed(siteBed, formatSiteBed.name, datasetId, clipSiteType)
    ## for peak bed, filter peak with  minPeakLen <= length <= maxPeakLen, and rename the peak_id
    ## combine formatted site bed and peak bed
    combineBed = os.path.join(outputDir, mergeBed)
    if siteFlag is True and peakFlag is True:
        CombineSitePeakBed(formatSiteBed.name, formatPeakBed.name, combineBed)
    elif siteFlag is False and peakFlag is True:
        CombineSitePeakBed('none', formatPeakBed.name, combineBed)
    elif siteFlag is True and peakFlag is False:
        CombineSitePeakBed(formatSiteBed.name, 'none', combineBed)
    ## delete temporary file
    if siteFlag is True:
        formatSiteBed.close()
        if peakFlag is True:
            kepetBed.close()
    if peakFlag is True:
        formatPeakBed.close()

########## main program ##########
outputDir = os.path.realpath(args.output)

## get all bed
bedList = list(map(os.path.realpath, glob(os.path.join(args.bed, "**", "*.bed"), recursive=True)))
bedDict = defaultdict(dict)
for bed in bedList:
    basename = os.path.basename(bed)
    bedDict[basename] = bed

clipSiteTypeDict = defaultdict(dict)
clipSiteTypeDict['PAR-CLIP'] = ['Mutation', 'peak']
clipSiteTypeDict['iCLIP'] = ['Truncation', 'peak']
clipSiteTypeDict['HITS-CLIP'] = ['Deletion', 'peak']
clipSiteTypeDict['eCLIP'] = ['Truncation', 'peak']
clipSiteTypeDict['irCLIP'] = ['Truncation', 'peak']
clipSiteTypeDict['pCLIP'] = ['Mutation', 'peak']
clipSiteTypeDict['PAR-iCLIP'] = ['Mutation', 'peak']
clipSiteTypeDict['FAST-iCLIP'] = ['Truncation', 'peak']
clipSiteTypeDict['CLEAR-CLIP'] = ['peak']
clipSiteTypeDict['PARCLIP-meRIP'] = ['Mutation', 'peak']
clipSiteTypeDict['sCLIP'] = ['peak']

clipSiteTypeShortDict = {'peak':'P', 'Mutation':'M', 'Truncation':'T', 'Deletion':'D'}
## record bed information from metadata table sheet
bedInfoDict = defaultdict(dict)
with open(args.meta, 'r') as f:
    __ = f.readline()
    for line in f:
        clipDict = {}
        row = line.strip().split('\t')
        datasetId = row[0]
        clipType = row[2]
        clipSource = row[6]
        mutation = row[-1]
        mergeBed = row[9]
        bedPrefix = mergeBed.replace('.bed', '')
        clipDict['peak'] = bedPrefix + '_peak.bed'
        clipDict['Mutation'] = bedPrefix + '_Mutation.bed'
        clipDict['Truncation'] = bedPrefix + '_Truncation.bed'
        clipDict['Deletion'] = bedPrefix + '_Deletion.bed'
        if mutation != 'none':
            clipDict['Mutation'] = bedPrefix + '_Mutation' + mutation + '.bed'
        clipSiteTypeList = clipSiteTypeDict[clipType]
        ##
        inforDict = defaultdict(dict)
        inforDict['datasetId'] = datasetId
        for clipSiteType in clipSiteTypeList:
            bedFileName = clipDict[clipSiteType]
            if bedFileName in bedDict:
                if clipSiteType == 'peak':
                    inforDict['peak'] = [clipSiteType, bedDict[bedFileName]]
                else:
                    inforDict['site'] = [clipSiteType, bedDict[bedFileName]]
        bedInfoDict[mergeBed] = inforDict

pool = Pool(processes=args.cpu)
with open(args.meta, 'r') as f:
    for mergeBed in sorted(bedInfoDict.keys()):
        #PoolBed(bedInfoDict, mergeBed, args.min, args.max, outputDir)
        pool.apply_async(PoolBed, args=(bedInfoDict, mergeBed, args.min, args.max, outputDir))
pool.close()
pool.join()
