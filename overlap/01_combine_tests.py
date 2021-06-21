#!/usr/bin/python -tt

import re
import sys
import os
from scipy import mean
from scipy import median
from scipy import stats
import numpy as np
import pybedtools as pbt
from collections import defaultdict 

def mapScaffsToChroms(scaffs2ChromFile):
    Chrom2Scaffs = defaultdict(list) # dict lists that maps chromosomes to scaffolds
    Scaffs2Chrom = defaultdict(str)
    scaffOrientation = {}
    scaffsInChroms = [] # list containing all scaffolds that map to a chromosome
    f = open(scaffs2ChromFile, 'r')
    for line in f:
        if not line.startswith("scaffoldName"):
            line = line.strip()
            info = line.split()
            chrom = info[1]
            scaff = info[0]
            scaffsInChroms.append(scaff)
            orientation = info[2]
            
            Chrom2Scaffs[chrom].append(scaff)
            Scaffs2Chrom[scaff] = chrom

            if orientation == "-":
                scaffOrientation[scaff] = "-"
            else:
                scaffOrientation[scaff] = "+"
    f.close()
    return(Chrom2Scaffs, Scaffs2Chrom, scaffOrientation, scaffsInChroms)

def getGeneBoundaries(geneBoundariesFile):
    geneBoundaries = {}
    f = open(geneBoundariesFile, 'r')
    for line in f:
        if not line.startswith("GeneName"):
            line = line.strip()
            info = line.split()
            geneName = info[0]
            geneNameEns = info[1]
            scaff = info[2]
            start = info[3]
            end = info[4]
            if scaff in geneBoundaries:
                geneBoundaries[scaff].append( (int(start), int(end), geneName, geneNameEns) )
            else:
                geneBoundaries[scaff]= []
                geneBoundaries[scaff].append( (int(start), int(end), geneName, geneNameEns) )
                    
    f.close()
    return(geneBoundaries)

def getSFStestsData(SFStests_file, scaffsAll, Chrom2Scaffs):

    fayWuHDict = defaultdict(list)
    TajDDict = defaultdict(list)

    fh = open(SFStests_file, 'r')
    for line in fh:
        if "window" not in line: 
            line = line.strip()
            line = line.split()
            scaff = line[0]
            winStart = int(line[2])-1
            winEnd = int(line[3])-1

            H = float(line[5])
            D = float(line[6])

            fayWuHDict[scaff].append( (winStart, winEnd, H) )
            TajDDict[scaff].append( (winStart, winEnd, D) )

    fh.close()
    return(fayWuHDict, TajDDict)


def getFstDataFromAngsd(fst_dir, scaffsAll, Chrom2Scaffs):
    fstDict = defaultdict(list)
    fstObs = 0
    for scaff in scaffsAll:
        fn = fst_dir + "/" + scaff + ".fst10kbwin2kbslide.txt" 
        if os.path.isfile(fn):
            f = open(fn, 'r')
            for line in f:
                if not line.startswith("region"):
                    line = line.strip()
                    line = line.split()
                    scaff = line[1]
                    fst = float(line[4])
                    sites = int(line[3])
                    if sites != 0:
                        # get window boundaries
                        winfo = line[0].split(")(")
                        winfo2 = winfo[2].split(",")
                        winfo2[1] = winfo2[1].replace(")","")
                        winStart = int(winfo2[0])-1
                        winEnd = int(winfo2[1])-1
                        fstDict[scaff].append( (winStart, winEnd, fst) )
                        fstObs += 1
                                    
            f.close()
        else:
            print("COULD NOT OPEN ANGSD FILE ", fn)
    return(fstDict, fstObs)	
	
def getDxyDataFromAngsd(dxy_dir, scaffsAll, Chrom2Scaffs):
    dxyDict = defaultdict(list)
    dxyObs = 0
    for scaff in scaffsAll:
        #fn = dxy_dir + "/Dxy_" + scaff + "_winSize10000.txt" 
        fn = dxy_dir + "/Dxy_" + scaff + "_winSize50_stepSize10.txt" 
        if os.path.isfile(fn):
            f = open(fn, 'r')
            for line in f:
                if not line.startswith("start"):
                    line = line.strip().split()
                    # get window boundaries
                    if line[0] != "NA" and line[2] != "NA":
                        winStart = int(line[0])-1
                        winEnd = int(line[1])-1
                        dxy = float(line[2])
                        dxyDict[scaff].append( (winStart, winEnd, dxy) )
                        dxyObs += 1
                                    
            f.close()
        else:
                print("COULD NOT OPEN ANGSD FILE ", fn)
    return(dxyDict, dxyObs)	

def getRAiSDdata(RAiSD_dir, scaffsAll, Chrom2Scaffs):
    muVarDict = defaultdict(list)
    muSfsDict = defaultdict(list)
    muLdDict = defaultdict(list)
    muCombDict = defaultdict(list)
    raisdObs = 0
    for scaff in scaffsAll:
        fn = RAiSD_dir + "/RAiSD_Report." + scaff 
        if os.path.isfile(fn):
            f = open(fn, 'r')
            for line in f:
                if scaff not in line:
                    line = line.strip()
                    line = line.split()
                    # get window boundaries subtract 1 for BED format
                    winStart = int(line[1]) - 1
                    winEnd = int(line[2]) - 1

                    var = float(line[3])
                    sfs = float(line[4])
                    ld = float(line[5])
                    comb = float(line[6])
                    
                    muVarDict[scaff].append( (winStart, winEnd, var) )
                    muSfsDict[scaff].append( (winStart, winEnd, sfs) )
                    muLdDict[scaff].append( (winStart, winEnd, ld) )
                    muCombDict[scaff].append( (winStart, winEnd, comb) )
                    raisdObs += 1
            f.close()
        else:
            print("COULD NOT OPEN RAiSD FILE ", fn)
    return(muVarDict, muSfsDict, muLdDict, muCombDict, raisdObs)	
    
def getHweData(HWE_file, minSegSites):
    hetDict = defaultdict(list)
    hetObs = 0
    f = open(HWE_file, 'r')
    for line in f:
        # skip header
        if "window" not in line:
            line = line.strip()
            line = line.split()
            # get window boundaries
            scaff = line[0]
            if line[2] != "NA":
                winStart = int(line[2]) - 1
                winEnd = int(line[3]) - 1
                hetDev = float(line[5])
                hetDict[scaff].append( (winStart, winEnd, hetDev) )
                hetObs += 1
    f.close()

    return(hetDict, hetObs)	

def getGene(scaff, winStart, winEnd, geneBoundaries):
    #geneBoundaries[scaff][geneName] = [ (int(start), int(end), geneName) ]
    geneName = "NA"
    geneNameEns = "NA"
    maxOverlap = 0
    maxOverlap_index = "NA"
    # some scaffolds, e.g. smaller ones, may not have genes!
    if scaff in geneBoundaries:
        for i in range(len(geneBoundaries[scaff])):
            overlap_bp = overlaps( (winStart,winEnd), (geneBoundaries[scaff][i][0],geneBoundaries[scaff][i][1]) )
            if overlap_bp > maxOverlap:
                maxOverlap_index = i
    if maxOverlap_index != "NA":
        geneName = geneBoundaries[scaff][maxOverlap_index][2]
        geneNameEns = geneBoundaries[scaff][maxOverlap_index][3]
    return(geneName, geneNameEns)

def reverseSorting(scaffOrientation, scaff):
    reverseTrueFalse = False
    if scaff in scaffOrientation:	
        if scaffOrientation[scaff] == "-":
            reverseTrueFalse = True
    return(reverseTrueFalse)

def makeArray(Dict, Obs):
    arr = np.empty(Obs)
    index = -1
    for scaff in Dict:
        for win in Dict[scaff]:
            index += 1
            winStart = win[0]
            winEnd = win[1]
            stat = win[2]
            arr[index] = stat
    return(arr)

def printOutlierList(name_prefix, Dict, cutoff, Scaffs2Chrom, percentileCutoff, outlier_dir):
    out_outlier = open(f"{outlier_dir}/{name_prefix}_{percentileCutoff}p_outliers.bed", 'w')
    out_all = open(f"{name_prefix}_all.bed", 'w')
    for scaff in Dict:
        chrom = Scaffs2Chrom[scaff]
        for win in Dict[scaff]:
            print(scaff, win[0], win[1], win[2], sep="\t", end="\n", file=out_all)
            if win[2] >= cutoff:
                print(scaff, win[0], win[1], win[2], sep="\t", end="\n", file=out_outlier)
    out_outlier.close()
    out_all.close()
def main():

    percentileCutoff = float(sys.argv[1])
    minSegSites = 29
    outlier_dir = "outliers"
    if not os.path.isdir(outlier_dir):
        os.mkdir(outlier_dir)

    ###########################################################################
    # DIRECTORIES AND FILES 
    ###########################################################################
    scaffs2ChromFile = "/tigress/bjarnold/AfricanElephants/genome/loxAfr4/ScaffoldsInChromosomes.txt"
    #geneBoundariesFile = "/n/holylfs/LABS/informatics/bjarnold/AfricanElephants/scripts/new/VCFanalyses/NewAnalysesWithOGs/position_2_gene.txt"

    # RAiSD dir
    RAiSD_dir = "/scratch/gpfs/bjarnold/AE/workflows/RAiSD/analysis"
    #RAiSD_dir = "/n/holylfs/LABS/informatics/bjarnold/AfricanElephants/scripts/new/VCFanalyses/NewAnalysesWithOGs/RAiSD/tusklessMZ/w50_c1"

    # FST FILES
    fst_dir = "/scratch/gpfs/bjarnold/AE/workflows/angsd/fst"

    # DXY FILES
    dxy_dir = "/scratch/gpfs/bjarnold/AE/workflows/angsd/dxy"

    ###
    # HET analyses
    ###
    # All tuskless alleles
    #HWE_file = "/scratch/gpfs/bjarnold/AE/workflows/HWE/HWEdeviations_Tuskless_25SNPwin_5stepSize_minAF0.15_maxAF0.85.txt"
    #HWE_file = "/scratch/gpfs/bjarnold/AE/workflows/HWE/HWEdeviations_Tuskless_50SNPwin_10stepSize_minAF0.1_maxAF0.9.txt"
    #HWE_file = "/scratch/gpfs/bjarnold/AE/workflows/HWE/HWEdeviations_Tuskless_50SNPwin_10stepSize_minAF0.0_maxAF1.0.txt"


    # Private tuskless alleles (i.e. not present in tusked samples)
    #HWE_file = "/tigress/bjarnold/AfricanElephants/data/HWE/TusklessPrivateAlleles/HWEdeviations_TusklessPrivateAlleles_50SNPwin_10stepSize_minAF0.1_maxAF0.9.txt"
    HWE_file = "/tigress/bjarnold/AfricanElephants/data/HWE/TusklessPrivateAlleles/HWEdeviations_TusklessPrivateAlleles_10SNPwin_2stepSize_minAF0.15_maxAF0.85.txt"
    #HWE_file = "/tigress/bjarnold/AfricanElephants/data/HWE/TusklessPrivateAlleles/HWEdeviations_TusklessPrivateAlleles_10SNPwin_2stepSize_minAF0.1_maxAF0.9.txt"
    #HWE_file = "/tigress/bjarnold/AfricanElephants/data/HWE/TusklessPrivateAlleles/HWEdeviations_TusklessPrivateAlleles_10SNPwin_2stepSize_minAF0.0_maxAF1.0.txt"

    ###########################################################################
    # GET INFO ABOUT GENES, SCAFFOLDS, CHROMOSOMES 
    ###########################################################################
    # FIND ALL SCAFFOLD NAMES 
    scaffsAll = []
    ih = open("/tigress/bjarnold/AfricanElephants/genome/scaffolds.txt",'r')
    for line in ih:
        line = line.split()
        s = line[0].strip(">")
        scaffsAll.append(s)
    ih.close()

    # FIND WHICH SCAFFOLDS MAP TO WHICH CHROMOSOMES
    #Chrom2Scaffs[chrom] = [scaffolds]
    #scaffOrientation[scaff] = "-"
    #scaffsInChroms = [scaffolds]
    Chrom2Scaffs, Scaffs2Chrom, scaffOrientation, scaffsInChroms = mapScaffsToChroms(scaffs2ChromFile)
    
    # FIND WHICH SCAFFOLDS DO NOT MAP TO CHROMOSOMES
    Chrom2Scaffs["Strays"] = []
    for scaff in scaffsAll:
        if scaff not in scaffsInChroms:
            Chrom2Scaffs["Strays"].append(scaff)

    # FIND GENE BOUNDARIES
    # geneBoundaries[scaff] = [ (int(start), int(end), geneNae) ]
    #geneBoundaries = getGeneBoundaries(geneBoundariesFile)

    ###########################################################################
    # LOAD RESULTS FOR OVERLAPPING
    ###########################################################################
    # FST RESULTS

    # raisdDict[scaff] = ( (winStart, winEnd, stat) )
    #raisdDict, raisdObs = getRAiSDdata(RAiSD_dir, scaffsAll, Chrom2Scaffs)
    muVarDict, muSfsDict, muLdDict, muCombDict, raisdObs = getRAiSDdata(RAiSD_dir, scaffsAll, Chrom2Scaffs)

    #fstDict[scaff] = ( (winStart, winEnd, fst) )   
    fstDict, fstObs = getFstDataFromAngsd(fst_dir, scaffsAll, Chrom2Scaffs)
    dxyDict, dxyObs = getDxyDataFromAngsd(dxy_dir, scaffsAll, Chrom2Scaffs)
    
    # hetDict[scaff] = ( (winStart, winEnd, hetDevTuskless) ) 
    hetDict, hetObs = getHweData(HWE_file, minSegSites)

    #fayWuHDict, TajDDict = getSFStestsData(SFStests_file, scaffsAll, Chrom2Scaffs)	

    muvar_arr = makeArray(muVarDict, raisdObs)
    musfs_arr = makeArray(muSfsDict, raisdObs)
    muld_arr = makeArray(muLdDict, raisdObs)
    mucomb_arr = makeArray(muCombDict, raisdObs)

    fst_arr = makeArray(fstDict, fstObs)
    dxy_arr = makeArray(dxyDict, dxyObs)
    het_arr = makeArray(hetDict, hetObs)

    print(np.percentile(mucomb_arr, percentileCutoff))
    print(np.percentile(fst_arr, percentileCutoff))
    print(np.percentile(dxy_arr, percentileCutoff))
    print(np.percentile(het_arr, percentileCutoff))

    muVar_cutoff = np.percentile(muvar_arr, percentileCutoff)
    muSfs_cutoff = np.percentile(musfs_arr, percentileCutoff)
    muLd_cutoff = np.percentile(muld_arr, percentileCutoff)
    muComb_cutoff = np.percentile(mucomb_arr, percentileCutoff)
    fst_cutoff = np.percentile(fst_arr, percentileCutoff)
    dxy_cutoff = np.percentile(dxy_arr, percentileCutoff)
    het_cutoff = np.percentile(het_arr, percentileCutoff)

    printOutlierList("RAiSD_muComb", muCombDict, muComb_cutoff, Scaffs2Chrom, percentileCutoff, outlier_dir)
    printOutlierList("RAiSD_muVar", muVarDict, muVar_cutoff, Scaffs2Chrom, percentileCutoff, outlier_dir)
    printOutlierList("RAiSD_muSfs", muSfsDict, muSfs_cutoff, Scaffs2Chrom, percentileCutoff, outlier_dir)
    printOutlierList("RAiSD_muLd", muLdDict, muLd_cutoff, Scaffs2Chrom, percentileCutoff, outlier_dir)
    printOutlierList("Fst", fstDict, fst_cutoff, Scaffs2Chrom, percentileCutoff, outlier_dir)
    printOutlierList("Dxy", dxyDict, dxy_cutoff, Scaffs2Chrom, percentileCutoff, outlier_dir)
    printOutlierList("Het", hetDict, het_cutoff, Scaffs2Chrom, percentileCutoff, outlier_dir)

if __name__ == '__main__':
  main()
