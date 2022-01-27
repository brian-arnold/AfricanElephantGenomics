#!/usr/bin/python -tt

import re
import sys
import os
import gzip
from statistics import mean 
from statistics import median
from collections import defaultdict

def parseMAFfile(fileName):
    afs = defaultdict(float)
    fh = gzip.open(fileName, 'rt', encoding='utf-8')
    for line in fh:
        if not line.startswith("chromo"):
            line = line.strip()
            line = line.split()
            pos = int(line[1])
            af = float(line[5])
            afs[pos] = af
    fh.close()
    return(afs)

def parsePOSfile(fileName, winSize, depthThresh):                                                                                                     
    posList = []                                                                                                                                  
    fh = gzip.open(fileName, 'rt', encoding='utf-8')                                                                                              
    for line in fh:                                                                                                                               
        if not line.startswith("chr"):                                                                                                        
            line = line.strip()                                                                                                           
            line = line.split()                                                                                                           
            pos = int(line[1])                                                                                                            
            depth = int(line[2])                                                                                                          
            if depth >= depthThresh:                                                                                                      
                posList.append(pos)
    fh.close()                                                                                                                                       
    return(posList)       

def calcDxy(list1, list2):
    dxyList = []
    if len(list1) != len(list2):
        return("NA")
        sys.exit("maf lists not equal!")
    for i in range(len(list1)):
        p1 = float(list1[i])	
        p2 = float(list2[i])	
        dxy = ((p1*(1-p2)) + (p2*(1-p1)))
        dxyList.append(dxy)
    # return the average
    return(dxyList)

def window(iterable, size, offset):
    # taken from stackoverflow and tested, except I added the offset and the try/except
    i = iter(iterable)
    win = []
    for e in range(0, size):
        try:
            win.append(next(i))
        except StopIteration:
            return
    yield win
    j = 0
    for e in i:
        win = win[1:] + [e]
        j += 1
        if j % offset == 0:
            yield win

def makeWindowedDict(List, winSize, offset):
    newDict = defaultdict(float)
    # List should be a (long) list that we will now dice up into windows
    # window returns an iterable
    wins = window(List, winSize, offset) 
    cnt = 1 # label windows numerically
    for win in wins:
        newDict[cnt] = mean(win)
        cnt += 1
    # were there enough SNPs for even a single window?
    if 1 not in newDict:
        newDict[1] = "NA"
    return(newDict)

def main():

    posFN = sys.argv[1] # list of all positions considered by ANGSD
    tusklessFN = sys.argv[2] # SNP sites in tuskless
    tuskedFN = sys.argv[3] # SNP sites in tusked

    winSize = 50 
    stepSize = 10 
    depthThresh = 56 # 14 individ * 4x depth = 56, THIS IS SAME NUMBER AS USED WHEN CALLING MAFs!!!
    scaffold = tusklessFN
    scaffold = scaffold.replace("_tuskless.mafs.gz", "")

    # first collect data into dicts, and use these to find positions that overlap in BOTH tuskless and tusked
    # also use sets for faster lookup?
    tusklessAFsDict = parseMAFfile(tusklessFN)
    tuskedAFsDict = parseMAFfile(tuskedFN)

    # find overlapping positions, and put these into lists that will then be used to make sliding window dicts
    positions = []
    tusklessAFsList = []
    tuskedAFsList = []
    for pos in sorted(tusklessAFsDict):
        if pos in tuskedAFsDict:
            positions.append(pos)
            tusklessAFsList.append( tusklessAFsDict[pos] )
    for pos in sorted(tuskedAFsDict):
        if pos in tusklessAFsDict:
            tuskedAFsList.append( tuskedAFsDict[pos] )
            
    # calculate dxy for all positions, store in to list
    dxyList = calcDxy(tusklessAFsList, tuskedAFsList)
    # use this list to make a windowed dict with mean dxy
    dxyWin = makeWindowedDict(dxyList, winSize, stepSize)

    # use positions list to make a windowed dict of start/end points of windows
    positionsWin = defaultdict(tuple)
    winSizes = []
    wins = window(positions, winSize, stepSize) 
    cnt = 1 # label windows numerically
    for win in wins:
        positionsWin[cnt] = (win[0], win[-1])
        winSizes.append(win[-1] - win[0])
        cnt += 1
    # were there enough SNPs for even a single window?
    if 1 not in positionsWin:
        positionsWin[1] = ("NA", "NA")

    if len(winSizes) > 0:
        print("mean win sizes: ", mean(winSizes))
        print("median win sizes: ", median(winSizes))

    oh = open(f"Dxy_{scaffold}_winSize{winSize}_stepSize{stepSize}.txt", 'w')
    print("start", "end", "dxy", file=oh)
    for win in sorted(dxyWin):
        start = positionsWin[win][0]
        end = positionsWin[win][1]
        print(start, end, dxyWin[win], file=oh)
    oh.close()

    sys.exit()

if __name__ == '__main__':
  main()
