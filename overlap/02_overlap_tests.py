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


def main():

    if len(sys.argv) < 2:
        sys.exit("specify a percentile")
    percentileCutoff = sys.argv[1]

    outlier_dir = "outliers"
    overlap_dir = "overlaps"
    if not os.path.isdir(overlap_dir):
        os.mkdir(overlap_dir)

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
    fstDxy = fst.intersect(dxy, wo=True).saveas(f"{overlap_dir}/FstDxy_overlap_{percentileCutoff}p.bed")

    # overlap muComb and muLd with FstDxy
    # wo = write original A and B entries plus # bp overlap between two features
    muCombFstDxy = muComb.intersect(fstDxy, wo=True).saveas(f"{overlap_dir}/muComb_FstDxy_overlap_{percentileCutoff}p.bed")
    muLdFstDxy = muLd.intersect(fstDxy, wo=True).saveas(f"{overlap_dir}/muLd_FstDxy_overlap_{percentileCutoff}p.bed")

    # overlap with het
    #het.intersect(muCombFstDxy, wo=True).saveas(f"{overlap_dir}/muComb_FstDxy_Het_overlap_{percentileCutoff}p.bed")
    hetMuLdFstDxy = het.intersect(muLdFstDxy, wo=True).saveas(f"{overlap_dir}/muLd_FstDxy_Het_overlap_{percentileCutoff}p.bed")

    muLdFstDxy.intersect(annot, wo=True).saveas(f"{overlap_dir}/muLd_FstDxy_overlap_{percentileCutoff}p_ANNOTATED.bed")
    hetMuLdFstDxy.intersect(annot, wo=True).saveas(f"{overlap_dir}/muLd_FstDxy_Het_overlap_{percentileCutoff}p_ANNOTATED.bed")
    muCombFstDxy.intersect(annot, wo=True).saveas(f"{overlap_dir}/muComb_FstDxy_overlap_{percentileCutoff}p_ANNOTATED.bed")

if __name__ == '__main__':
  main()
