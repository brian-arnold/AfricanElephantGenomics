#!/usr/bin/python -tt

import re
import sys
import os
import gzip
from collections import defaultdict

def main():

    chromToAnalyze = "1"
    scaffs2ChromFile = "/tigress/bjarnold/AfricanElephants/genome/loxAfr4/ScaffoldsInChromosomes.txt"
    maxDist = 200000
    r2ByDist = {}
    r2ByDist_totals = {}
    # initialize
    for i in range(maxDist):
        r2ByDist[i+1] = 0.0
        r2ByDist_totals[i+1] = 0
    
    scaffsInChroms = []
    f = open(scaffs2ChromFile,'r')
    for line in f:
        if not line.startswith("scaffoldName"):
            line = line.strip()
            info = line.split()
            chrom = info[1]
            scaff = info[0]
            if chrom == chromToAnalyze:
                scaffsInChroms.append(scaff)
    f.close()

    for scaffold in scaffsInChroms:
        fileName = scaffold + ".geno.ld"
        f = open(fileName, 'r')
        for line in f:
            #line = line.decode("utf-8") # need to convert bytes to string when working with compressed files
            line = line.strip()
            data = line.split('\t')
            if data[0] != 'CHR' and data[4] != '-nan': 
                dist = int(data[2]) - int(data[1])
                if dist <= maxDist:
                    r2ByDist[dist] += float(data[4])
                    r2ByDist_totals[dist] += 1
        f.close()
    outName = "Chr" + chromToAnalyze + "_r2ByDist.txt"
    o = open(outName,'w')
    for dist in range(maxDist):
        dist=dist+1
        if r2ByDist_totals[dist] > 0:
            print(dist, "\t", r2ByDist_totals[dist], "\t", r2ByDist[dist]/r2ByDist_totals[dist], file=o)
    o.close()

if __name__ == '__main__':
  main()
