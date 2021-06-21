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
    Scaffs2Chrom = defaultdict(list) # dict lists that maps chromosomes to scaffolds
    f = open(scaffs2ChromFile, 'r')                                                                                                               
    for line in f:                                                                                                                                
        if not line.startswith("scaffoldName"):                                                                                               
            line = line.strip()                                                                                                           
            info = line.split()                                                                                                           
            chrom = info[1]                                                                                                               
            scaff = info[0]                                                                                                               
            orientation = info[2]                                                                                                         
            Chrom2Scaffs[chrom].append(scaff)                                                                                     
            Scaffs2Chrom[scaff] = chrom
    f.close()                                                                                                                                     
    return(Chrom2Scaffs, Scaffs2Chrom)      

def getDentalGenes(fn):
    f = open(fn, 'r')
    dentalGenes = []
    for line in f:
        line = line.strip().split("\t")
        genes = line[1].split(",")
        for i in genes:
            if i != "NA":
                dentalGenes.append(i)
    return(dentalGenes)

def main():

    percentileCutoff = sys.argv[1]
    scaffs2ChromFile = "/tigress/bjarnold/AfricanElephants/genome/loxAfr4/ScaffoldsInChromosomes.txt"  
    Chrom2Scaffs, Scaffs2Chrom = mapScaffsToChroms(scaffs2ChromFile)

    dentalGenes = getDentalGenes("/tigress/bjarnold/AfricanElephants/data/DentalGenesWithEnsLafgIDs.txt")

    ###
    # OVERALLPING HET, MULD, FST, DXY
    ###
    overlapFile = f"overlaps/muLd_FstDxy_Het_overlap_{percentileCutoff}p_ANNOTATED.bed"

    o = open(f"muLd_FstDxy_Het_overlap_{percentileCutoff}p_toothGenes.bed", 'w')
    print("chrom", "scaff", "winStart", "winEnd", "dxy", "fst", "muLd", "het", "gene", "geneEnsembl", file=o)
    f = open(overlapFile, 'r')
    for line in f:
        line = line.strip()
        info = line.split()
        scaff = info[0]
        winStart = [ info[1], info[5], info[9], info[13] ] 
        winEnd = [ info[2], info[6], info[10], info[14] ] 
        dxy = info[15]
        fst = info[11]
        muLd = info[7]
        het = info[3]
        gene = info[22]
        geneEnsembl = info[23]

        p = 0
        for i in dentalGenes:
            if i in line:
                p=1
        if p:
            winStart = [int(i) for i in winStart]
            winEnd = [int(i) for i in winEnd]
            print(Scaffs2Chrom[scaff], scaff, min(winStart), max(winEnd), dxy, fst, muLd, het, gene, geneEnsembl, file=o)
    o.close()
    f.close()

if __name__ == '__main__':
  main()
