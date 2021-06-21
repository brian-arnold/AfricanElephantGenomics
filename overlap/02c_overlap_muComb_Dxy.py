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


    overlap_dir = "overlaps"
    if not os.path.isdir(overlap_dir):
        os.mkdir(overlap_dir)

    muComb_file = f"RAiSD_muComb_all.bed"
    dxy_file = f"Dxy_all.bed"

    annot_file = "/tigress/bjarnold/AfricanElephants/scripts/annotations/position_2_gene.bed"

    muComb = pbt.BedTool(muComb_file)
    dxy = pbt.BedTool(dxy_file)
    annot = pbt.BedTool(annot_file)

    # wo = write original A and B entries plus # bp overlap between two features
    muCombDxy = muComb.intersect(dxy, wo=True).saveas(f"{overlap_dir}/muComb_Dxy_overlap_All.bed")

if __name__ == '__main__':
  main()
