#!/usr/bin/python -tt

import re
import sys
import os

def main():
    """
    This script uses the Chromosomes.v2.map from Palkopoulou et al (2018), PNAS and converts it into a table
    that lists which scaffolds belong to which chromosome along with the orientation of each scaffold on 
    each chromosome.
    Admittedly, I had to make an educated guess as to what various columns of the Chromosomes.v2.map actually
    represented, as the authors of the study that produced this file did not respond to emails. Nonetheless,
    information in this file was consistent with known scaffold lengths, based on the 3.0 Loxodonta africana
    genome assembly.
    """

    f = open("Chromosomes.v2.map", 'r')
    o = open("ScaffoldsInChromosomes.txt", 'w')
    print("scaffoldName\tchromosome\torientation", file=o)
    ChromScaff = {}
    MarkersToSkip = {}
    contigList = {}
    for line in f:
        line = line.strip()
        info = line.split()
        # FILTER OUT BAD MARKERS
        if line.startswith("ManualTag"):
            MarkersToSkip[info[3]] = 1
        if line.startswith("chr"):
            if info[0] == "chr" and info[1] == "marker":
                count = 0
                for line2 in f:
                    if line2.startswith("Synopsis"):
                        break
                    else:
                        count += 1
                        line2 = line2.strip()
                        info2 = line2.split()
                        chrom = info2[0]
                        contig = int(info2[3])
                        scaff = info2[4]
                        markerID = info2[1]
                        if markerID not in MarkersToSkip:
                            ChromScaff[scaff] = (chrom, count)
                            if scaff in contigList:
                                contigList[scaff].append(contig)
                            else:
                                contigList[scaff] = []
                                contigList[scaff].append(contig)
    f.close()
    #for scaff in ChromScaff:
    for scaff in sorted(ChromScaff, key=lambda k: (str(ChromScaff[k][0]), int(ChromScaff[k][1])), reverse=False):
            scaffName = "scaffold_" + str(scaff)
            print(scaffName, ChromScaff[scaff][0], sep="\t", end="\t", file=o)
            summ = 0
            for i in range(len(contigList[scaff])-1):
                summ += contigList[scaff][i+1] - contigList[scaff][i] 
            if summ > 0:
                print("+", file=o)
            elif summ < 0:
                print("-", file=o)
            else:
                print("NA", file=o)

if __name__ == '__main__':
  main()
