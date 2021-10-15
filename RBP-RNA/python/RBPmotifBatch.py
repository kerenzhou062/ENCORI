#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
import datetime
import shutil
from glob import glob
import subprocess
import multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument('-cpu', action='store', type=int,
                    default=10, help='cores for multiprocessing')
parser.add_argument('-fasta', action='store', type=str,
                    help='The genome Fasta')
parser.add_argument('-anno', action='store', type=str,
                    help='The genome annotation bed12')
parser.add_argument('-genome', action='store', type=str,
                    help='The genome(hg19, mm10, ce10, dm6)')
parser.add_argument('-input', action='store', type=str,
                    help='The input bed folder')
parser.add_argument('-ID', nargs='*', type=str,
                    help='The designated sampleid list')
parser.add_argument('-meta', action='store', type=str,
                    help='metadata matrix of CLIP-seq')
parser.add_argument('-agoRef', action='store', type=str,
                    help='The AGO reference file')
parser.add_argument('-length', nargs='+', type=int,
                    default=[4, 6, 8],
                    help='the length used for homer de novo motif finding')
parser.add_argument('-motifRef', action='store', type=str,
                    help='The motif reference from ATtRACT_db')
parser.add_argument('-size', action='store', type=str,
                    default='200', help='The -size of findMotifsGenome.pl')
parser.add_argument('-output', action='store', type=str,
                    default="cluster.txt", help='The output cluster file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

def RunMotif(command, datasetID, bedInput, outputDir, preparseDir, logFile):
    # delete the existing target folder
    if os.path.isdir(outputDir):
        shutil.rmtree(outputDir)
    os.makedirs(outputDir)
    ## change the working dir
    os.chdir(outputDir)

    bedDict = defaultdict(list)
    with open(bedInput, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            chro = row[0]
            center = int((int(row[1]) + int(row[2])) / 2)
            posID = row[3]
            strand = row[5]
            bedDict[posID].extend([chro, center, strand])
    stdout = open(logFile, 'w')
    subprocess.call(command, stdin=None, stdout=stdout, stderr=stdout, shell=True)
    stdout.close()
    shutil.rmtree(preparseDir)
    motifMatrixList = glob(os.path.join('homerResults', '*.motif'))
    motifBedDir = 'motifBed'
    os.makedirs(motifBedDir)
    motifBedListFile = open('motifBedList.txt', 'w')
    for motifMatrix in sorted(motifMatrixList):
        matrixName = os.path.basename(motifMatrix)
        if matrixName in bedMotifLib:
            basename = os.path.splitext(matrixName)[0]
            # find motif in bed and transform to temp-bed format
            findMotifCommand = 'findMotifsGenome.pl {0} {1} {2} -find {3} -size {4} -rna'
            findMotifCommand = findMotifCommand.format(bedInput, args.genome, '.', motifMatrix, args.size)
            motifTemp1 = os.path.join(motifBedDir, datasetID + '-' + basename + '.tmp1')
            motifTemp2 = os.path.join(motifBedDir, datasetID + '-' + basename + '.tmp2')
            motifBed = os.path.join(motifBedDir, datasetID + '-' + basename + '.bed')
            motifResult = subprocess.check_output(findMotifCommand, stdin=None, stderr=subprocess.DEVNULL, shell=True)
            motifResultList = bytes.decode(motifResult).split('\n')[1:]
            motifResultDict = defaultdict(list)
            for line in motifResultList:
                row = line.strip().split('\t')
                if len(row) != 6:
                    continue
                posID = row[0]
                offset = int(row[1])
                motifSeq = row[2]
                motifScore = row[-1]
                bedInfo = bedDict[posID]
                chro = bedInfo[0]
                center = bedInfo[1]
                strand = bedInfo[2]
                if strand == '+':
                    start = center + offset - 1
                    end = start + len(motifSeq)
                else:
                    end = center - offset
                    start = end - len(motifSeq)
                key = ':'.join([chro, str(start), str(end), strand])
                motifBedList = [posID, motifSeq, motifScore]
                motifResultDict[key].append(motifBedList)
            # writing data to temp-bed
            with open(motifTemp1, 'w') as temp:
                for key in sorted(motifResultDict.keys()):
                    tempList = key.split(':')
                    motifSeq = motifResultDict[key][0][1]
                    tempBedList = tempList[0:3]
                    tempBedList.append(key + '|' + motifSeq)
                    tempBedList.extend(['0', tempList[3]])
                    temp.write('\t'.join(tempBedList) + '\n')
            # get fasta from temp-bed and check the motifseq
            getfastaComand = 'bedtools getfasta -fi {0} -bed {1} -s -tab | sort | uniq'.format(args.fasta, motifTemp1)
            getfastaResult = subprocess.check_output(getfastaComand, stdin=None, stderr=subprocess.DEVNULL, shell=True)
            getfastaResultList = bytes.decode(getfastaResult).split('\n')
            with open(motifTemp2, 'w') as temp:
                for line in getfastaResultList:
                    row = line.strip().split('\t')
                    if len(row) != 2:
                        continue
                    fastaIDList = re.split(r'\(|:', row[0])
                    locusList = fastaIDList[1].split('-')
                    fastaIDList[-1] = fastaIDList[-1].replace(')', '')
                    keyList = [fastaIDList[0], locusList[0], locusList[1], fastaIDList[-1]]
                    key = ':'.join(keyList)
                    fastaSeq = row[1].upper()
                    motifBedBulkyList = motifResultDict[key]
                    posIDlist = list()
                    motifScore = 0
                    for motifBedList in motifBedBulkyList:
                        motifSeq = motifBedList[1]
                        if motifSeq != fastaSeq:
                            continue
                        posIDlist.append(motifBedList[0])
                        if float(motifBedList[2]) > motifScore:
                            motifScore = float(motifBedList[2])
                    if posIDlist:
                        posIDcat = ':'.join(posIDlist)
                        motifID = '|'.join([posIDcat, fastaSeq])
                        bedLineList = keyList[0:3]
                        bedLineList.extend([motifID, str(motifScore), keyList[3]])
                        temp.write('\t'.join(bedLineList) + '\n')
            annoComand = 'bedtools intersect -a {0} -b {1} -wa -wb -s'.format(motifTemp2, args.anno)
            annoResult = subprocess.check_output(annoComand, stdin=None, stderr=subprocess.DEVNULL, shell=True)
            annoResultList = bytes.decode(annoResult).split('\n')
            motifBedDict = defaultdict(dict)
            for line in annoResultList:
                if line == '':
                    continue
                row = line.strip().split('\t')
                bedID = row[3]
                motifBedDict[bedID]['bed'] = row[0:6]
                geneInfo = '|'.join(row[9].split('|')[3:6])
                if 'gene' in motifBedDict[bedID]:
                    motifBedDict[bedID]['gene'].append(geneInfo)
                else:
                    motifBedDict[bedID]['gene'] = [geneInfo]

            # delete tempfile
            os.remove(motifTemp1)
            os.remove(motifTemp2)
            # write motif-bed
            with open(motifBed, 'w') as out:
                for bedID in sorted(motifBedDict.keys()):
                    bedRow = motifBedDict[bedID]['bed']
                    genes = ','.join(sorted(set(motifBedDict[bedID]['gene'])))
                    bedRow[3] = '='.join([bedRow[3], genes])
                    out.write('\t'.join(bedRow) + '\n')
            compressComand = 'gzip {0}'.format(motifBed)
            subprocess.run(compressComand, stdin=None, stderr=subprocess.DEVNULL, shell=True)
            motifBedListFile.write(os.path.basename(motifBed) + '\n')
    motifBedListFile.close()

    # delete tempFile
    tempFileList = [glob(os.path.join(outputDir, '*.tmp')), glob(os.path.join(outputDir, '*.pos'))]
    for tempFiles in tempFileList:
        if tempFiles:
            for tempFile in tempFiles:
                if os.path.isfile(tempFile):
                    os.remove(tempFile)


starttime = datetime.datetime.now()
sys.stderr.write("Program starting!\n")

bedMotifLib = ['motif' + str(i + 1) + '.motif' for i in range(25)]

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

metaDict = defaultdict(dict)
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
            metaDict[datasetID]['rbp'] = RBP.upper()
            metaDict[datasetID]['file'] = bedDict[fileName]

motifDict = defaultdict(list)
if bool(args.motifRef):
    with open(args.motifRef, 'r') as f:
        #geneName,motifSeq,length
        __ = f.readline()
        for line in f:
            row = line.strip().split('\t')
            motifDict[row[0]].append(int(row[2]))

pool = multiprocessing.Pool(processes=args.cpu)

if args.ID:
    runIDlist = sorted(args.ID)
else:
    runIDlist = sorted(metaDict.keys())

for datasetID in runIDlist:
    RBP = metaDict[datasetID]['rbp']
    bedInput = metaDict[datasetID]['file']
    if RBP in motifDict:
        refLenAve = round(sum(motifDict[RBP]) / len(motifDict[RBP]))
        if refLenAve < 5:
            refLenAve = 5
        elif refLenAve > 8:
            refLenAve = 8
        refLenArr = list(map(str, [refLenAve-2, refLenAve, refLenAve + 2]))
    else:
        refLenArr = list(map(str, args.length))
    # preparing log, output folder, preparse folder
    outputDir = os.path.join(args.output, datasetID)
    logFile = os.path.join(outputDir, datasetID + '.log')
    preparseDir = os.path.join(outputDir, 'preparsed')
    refLen = ','.join(refLenArr)
    command = 'findMotifsGenome.pl {0} {1} {2} -rna -noknown -size {3} -len {4} -p 1 -S 25 -preparse -preparsedDir {5}'
    command = command.format(bedInput, args.genome, outputDir, args.size, refLen, preparseDir)
    #RunMotif(command, datasetID, bedInput, outputDir, preparseDir, logFile)
    pool.apply_async(RunMotif, (command, datasetID, bedInput, outputDir, preparseDir, logFile))

pool.close()
pool.join()

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("All jobs done!")
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
