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


def make_venn( a, b, c, d, l  ):
    # find overlaps, print them in sorted order so that overlaps have same string across files
    str_vec = []
    for i in range( (a-b-c-d).count() ):
       str_vec.append(l[0]) 
    for i in range( (a+b-c-d).count() ):
       x = "_".join( sorted([ l[0], l[1], str(i) ]) )
       str_vec.append(x) 
    for i in range( (a-b+c-d).count() ):
       x = "_".join( sorted([ l[0], l[2], str(i) ]) )
       str_vec.append(x) 
    for i in range( (a-b-c+d).count() ):
       x = "_".join( sorted([ l[0], l[3], str(i) ]) )
       str_vec.append(x) 
    for i in range( (a+b+c-d).count() ): 
       x = "_".join( sorted([ l[0], l[1], l[2], str(i) ]) )
       str_vec.append(x) 
    for i in range( (a+b-c+d).count() ):
       x = "_".join( sorted([ l[0], l[1], l[3], str(i) ]) )
       str_vec.append(x) 
    for i in range( (a-b+c+d).count() ):
       x = "_".join( sorted([ l[0], l[2], l[3], str(i) ]) )
       str_vec.append(x) 
    for i in range( (a+b+c+d).count() ):
       x = "_".join( sorted([ l[0], l[1], l[2], l[3], str(i) ]) )
       str_vec.append(x) 
    return(str_vec)

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

    if len(sys.argv) < 2:
        sys.exit("specify a percentile")
    percentileCutoff = sys.argv[1]

    outlier_dir = "outliers"
    overlap_dir = "overlaps"
    if not os.path.isdir(overlap_dir):
        os.mkdir(overlap_dir)

    X_scaffs = ["scaffold_94", "scaffold_39", "scaffold_56", "scaffold_24", "scaffold_32", "scaffold_120"]
    dentalGenes = getDentalGenes("/tigress/bjarnold/AfricanElephants/data/DentalGenesWithEnsLafgIDs.txt")

    muComb_file = f"{outlier_dir}/RAiSD_muComb_{percentileCutoff}p_outliers.bed"
    muVar_file = f"{outlier_dir}/RAiSD_muVar_{percentileCutoff}p_outliers.bed"
    muSfs_file = f"{outlier_dir}/RAiSD_muSfs_{percentileCutoff}p_outliers.bed"
    muLd_file = f"{outlier_dir}/RAiSD_muLd_{percentileCutoff}p_outliers.bed"

    fst_file = f"{outlier_dir}/Fst_{percentileCutoff}p_outliers.bed"
    dxy_file = f"{outlier_dir}/Dxy_{percentileCutoff}p_outliers.bed"
    het_file = f"{outlier_dir}/Het_{percentileCutoff}p_outliers.bed"

    annot_file = "/tigress/bjarnold/AfricanElephants/scripts/annotations/position_2_gene.bed"

    muComb = pbt.BedTool(muComb_file)
    muVar = pbt.BedTool(muVar_file)
    muSfs = pbt.BedTool(muSfs_file)
    muLd = pbt.BedTool(muLd_file)

    fst = pbt.BedTool(fst_file)
    dxy = pbt.BedTool(dxy_file)
    het = pbt.BedTool(het_file)

    annot = pbt.BedTool(annot_file)

    # overlap Fst and Dxy
    # for testing and making sure make_venn() func works
    print(fst.count())
    print((fst-dxy-het-muLd).count())
    print((fst+dxy-het-muLd).count())
    print((fst+dxy+het+muLd).count())

    # Uncomment the following block to make lists for the 4-way venn diagram
    # with fst, dxy, LD, het
    fst_str = make_venn( fst, dxy, muLd, het, ["fst", "dxy", "muLD", "het"]  )
    o = open("fst_venn.txt", 'w')
    for i in fst_str:
        print(i, file=o)

    dxy_str = make_venn( dxy, fst, muLd, het, ["dxy", "fst", "muLD", "het"]  )
    o = open("dxy_venn.txt", 'w')
    for i in dxy_str:
        print(i, file=o)

    muLd_str = make_venn( muLd, fst, dxy, het, ["muLD", "fst", "dxy", "het"]  )
    o = open("muLd_venn.txt", 'w')
    for i in muLd_str:
        print(i, file=o)

    het_str = make_venn( het, fst, dxy, muLd, ["het", "fst", "dxy", "muLD"]  )
    o = open("het_venn.txt", 'w')
    for i in het_str:
        print(i, file=o)


    ##### THIS LOOKS AT ALL THE OVERLAPS FOR ALL 4 STATS, AND SEES WHICH ARE ON X, THEN WHICH ARE TEETH GENES
    all_intersect = (het+fst+dxy+muLd)
    all_intersect_annot = all_intersect.intersect(annot, wo=True)
    print("number of windows overlapping all 4 stats, sanity check:", len(all_intersect))

    # filter for only windows in X chromosome, save as file
    X_windows = 0
    X_all_overlaps = f"All_{percentileCutoff}p_outliers_X.bed"
    out = open(X_all_overlaps, 'w')
    for i in all_intersect:
        scaff = i[0]
        if scaff in X_scaffs:
            X_windows += 1
            print(i, end="", file=out)
    print("Num 4X overlap windows on X:", X_windows)

    # filter for windows on X, ANNOTATED
    X_all_overlaps = f"All_{percentileCutoff}p_outliers_X_ANNOTATED.bed"
    out = open(X_all_overlaps, 'w')
    for i in all_intersect_annot:
        scaff = i[0]
        if scaff in X_scaffs:
            print(i, end="", file=out)

    # filter for windows on X, ANNOTATED, TEETH
    X_all_overlaps = f"All_{percentileCutoff}p_outliers_X_ANNOTATED_toothGenes.bed"
    out = open(X_all_overlaps, 'w')
    for i in all_intersect_annot:
        scaff = i[0]
        if scaff in X_scaffs:
            for j in i:
                for k in dentalGenes:
                    if k in j:
                        print(i, end="", file=out)

    # get all annotated windows, not just those on X
    all_overlaps = f"All_{percentileCutoff}p_outliers_ANNOTATED.bed"
    out = open(all_overlaps, 'w')
    for i in all_intersect_annot:
        print(i, end="", file=out)

    sys.exit()

if __name__ == '__main__':
  main()
