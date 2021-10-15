#!/usr/bin/env python3
import os
import sys
import re
import argparse
from collections import defaultdict
import tempfile
import subprocess
import datetime
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-input', action='store', type=str,
                    help='The mirTarget file')
parser.add_argument('-conscore', action='store', type=str,
                    help='conservative score bedgraph file')
parser.add_argument('-output', action='store', type=str,
                    default="miRmRNA.txt", help='The output file')

args = parser.parse_args()
if len(sys.argv[1:])==0:
        parser.print_help()
        parser.exit()

starttime = datetime.datetime.now()

chrRegex = re.compile(r'^chr(\d+|[XYMxym])')
mirBedList = []
with open(args.input, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        if bool(chrRegex.match(row[1])) is True:
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            mirBedList.append([chrom, start, end])

mirBedList = sorted(mirBedList, key=lambda x:(x[0], x[1]))

mirBed6Tmp = tempfile.NamedTemporaryFile(suffix='.tmp', delete=True)
with open(mirBed6Tmp.name, 'w') as out:
    for mirBedRow in mirBedList:
        out.write('\t'.join(map(str, mirBedRow)) + '\n')

command = 'bedtools merge -i {} | sort -k1,1 -k2,2n | \
  bedtools intersect -a {} -b stdin -u -sorted -nonamecheck'.format(mirBed6Tmp.name, args.conscore)

filterResList = bytes.decode(subprocess.check_output(command, shell=True)).split('\n')
mirBed6Tmp.close()

consDict = defaultdict(dict)
for line in filterResList:
    if bool(line) is True:
        row = line.strip().split('\t')
        chrom = row[0]
        start = int(row[1])
        end = int(row[2])
        coorList = list(range(start, end))
        score = float(row[-1])
        for coor in coorList:
            consDict[chrom][coor] = score

with open(args.output, 'w') as out:
    with open(args.input, 'r') as f:
        ## header
        row = f.readline().strip().split('\t')
        row += ['consAve\n']
        out.write('\t'.join(row))
        for line in f:
            row = line.strip().split('\t')
            chrom = row[1]
            start = int(row[2])
            end = int(row[3])
            coorList = list(range(start, end))
            scoreList = []
            for coor in coorList:
                if coor in consDict[chrom]:
                    scoreList.append(consDict[chrom][coor])
                else:
                    scoreList.append(0)
            consAve =  np.average(scoreList)
            row += [str(round(consAve, 3)) + '\n']
            out.write('\t'.join(row))

endtime = datetime.datetime.now()
collapsed = (endtime - starttime).seconds
sys.stderr.write("Total collapsed time: {0}s\n".format(collapsed))
