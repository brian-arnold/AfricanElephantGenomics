#!/usr/bin/python -tt

import re
import sys
import os
from cyvcf2 import VCF
from statistics import mean
from collections import defaultdict
# see https://brentp.github.io/cyvcf2/docstrings.html

#def getRGinfo(fileName):

def filterVariantMetrics(variant, Filter):
    PASS=1
    for key in Filter:
        if variant.INFO.get(key) == None:
            PASS=0
    if PASS: # site hasn't failed yet, all filter criteria available
        if variant.INFO.get('QD') < Filter['QD']:
            PASS=0
        if variant.INFO.get('FS') > Filter['FS']:
            PASS=0
        if variant.INFO.get('SOR') > Filter['SOR']:
            PASS=0
        if variant.INFO.get('MQ') < Filter['MQ']:
            PASS=0
        if variant.INFO.get('MQRankSum') < Filter['MQRankSum']:
            PASS=0
        if variant.INFO.get('ReadPosRankSum') < Filter['ReadPosRankSum']:
            PASS=0
        #if variant.INFO.get('ExcessHet') > Filter['ExcessHet']:
        #	PASS=0
                #print(variant.gt_types)
    return(PASS)

def filterVariantDepthMissData(gt_depths, gt_types, Filter):
    PASS = 1
    missing = 0
    for i in range(len(gt_depths)):
        if (gt_depths[i] < Filter['MinDepthPerInd']) or (gt_depths[i] > Filter['MaxDepthPerInd']) or (gt_types[i] == 3): # has low depth or missing genotype (unlikely for reasonable depth thresholds), but for which category?
            missing += 1
    if missing > Filter['MaxMissIndPerCategory']: 
        PASS = 0
    return(PASS)

def CalcHWEdev(gt_types):
    n = len(gt_types)
    AF = 0.0
    #obs
    homRef_obs = 0.0
    het_obs = 0.0
    homAlt_obs = 0.0
    for i in range(n):	
        gt = gt_types[i]
        if gt == 0:
            homRef_obs += 1/n
        elif gt == 1:
            het_obs += 1/n
            AF += 0.5/n
        elif gt == 2:
            homAlt_obs += 1/n
            AF += 1/n
        else:
            print("oh noooo")

    homRef_exp = (1-AF)*(1-AF)
    het_exp = 2*AF*(1-AF)
    homAlt_exp = AF*AF

    HWEdev = abs(homRef_obs - homRef_exp) + abs(het_obs - het_exp) + abs(homAlt_obs - homAlt_exp) 
    HETdev = (het_obs - het_exp)
    HOMdev = (homRef_obs - homRef_exp) + (homAlt_obs - homAlt_exp)
    return(HWEdev, HETdev, HOMdev)

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

def makeWindowedDict(Dict, winSize, offset):
    newDict = defaultdict(lambda: defaultdict(float))
    for scaff in Dict:
        # Dict[scaff] should be a (long) list that we will now dice up into windows
        #if len(Dict[scaff]) >= winSize:         # if scaffold doesn't have enough SNPs for full window, throws runtime error
        wins = window(Dict[scaff], winSize, offset) 
        cnt = 1 # label windows numerically
        for win in wins:
            newDict[scaff][cnt] = mean(win)
            cnt += 1
        #else:
        #    newDict[scaff][1] = "NA"
        # if window returned empty iterable
        if scaff not in newDict:
            newDict[scaff][1] = "NA" 
    return(newDict)

def main():

    winSize = int(sys.argv[1]) # the size of the SNP window
    offset = int(sys.argv[2])
    minAF_tuskless = 0.15 # the minimum allele frequency in tuskless samples
    maxAF_tuskless = 0.85 # the maximum allele frequency in tuskless samples
    minAF_tusked = 0.00 # the minimum allele frequency in tusked samples
    maxAF_tusked= 0.00 # the maximum allele frequency in tusked samples

    vcf = VCF('/scratch/gpfs/bjarnold/AE/VCF/Final_filtered_annotated.vcf', gts012=True) # gts012=True makes value of 3 UNKNOWN, with 0,1,2 corresponding to numb ALT alleles
    #vcf = VCF('/scratch/gpfs/bjarnold/AE/VCF/test_100k.vcf', gts012=True) # gts012=True makes value of 3 UNKNOWN, with 0,1,2 corresponding to numb ALT alleles

    Groups = {} # 0 = tuskless, 1 = tusked
    Groups[0] = ["0045B", "2983B", "2984B", "G17A", "G19A", "G22A", "T2B"] # T1A excluded, contaminated
    Groups[1] = ["2981B", "2982B", "2985B", "G18A", "G20A", "G21A"] 
    Samples = Groups[0] + Groups[1]
    vcf.set_samples(Samples) 

    two_n = len(Samples)*2 # for diploid, for calculating statistics
    two_n_tuskless = len(Groups[0])*2 # for diploid, for calculating statistics
    two_n_tusked = len(Groups[1])*2 # for diploid, for calculating statistics

    #############
    # for each index in vcf.samples list, which group does it correspond to
    #############
    index2Group = {}
    for index in range(len(vcf.samples)):
        for grp in Groups:
            if vcf.samples[index] in Groups[grp]:
                index2Group[index] = grp   

    #############
    # construct dict with filter metrics
    #############
    FilterMetrics = {}
    FilterMetrics['QD'] = 2 # fail variants less than 2
    FilterMetrics['FS'] = 60 # fail variants greater than 60
    FilterMetrics['SOR'] = 3 # fail variants greater than 3
    FilterMetrics['MQ'] = 40 # fail variants less than 40
    FilterMetrics['MQRankSum'] = -12.5 # fail variants less than -12.5 
    FilterMetrics['ReadPosRankSum'] = -8 # fail variants less than -8 
    FilterDP = {}
    FilterDP['MinDepthPerInd'] = 8
    FilterDP['MaxDepthPerInd'] = 80
    FilterDP['MaxMissIndPerCategory'] = 0

    HWEdev = defaultdict(list)
    HETdev = defaultdict(list)
    HOMdev = defaultdict(list)
    depthsDict = defaultdict(list)
    positions = defaultdict(list) 

    currWin = 0
    snpCount = 0

    for variant in vcf:
        # only look at biAllelic Sites
        if len(variant.ALT) == 1 and variant.var_type == "snp":
            PASS_metrics = filterVariantMetrics(variant, FilterMetrics)
            # go thru gt_depths, check if greater than Filter, tally up miss data for each category
            PASS_depthMissData = filterVariantDepthMissData(variant.gt_depths, variant.gt_types, FilterDP)

            # gt_types is array of 0,1,2,3==HOM_REF, HET, HOM_ALT, UNKNOWN 
            if PASS_metrics == 1 and PASS_depthMissData == 1:
                # collect allele frequency info for further filtering
                AF = 0
                AF_tuskless = 0
                AF_tusked = 0
                # Get AF for tuskless and tusked groups
                for index in range(len(variant.gt_types)):
                    AF += variant.gt_types[index]
                    if index2Group[index] == 0:
                        AF_tuskless += variant.gt_types[index] 
                    else:
                        AF_tusked += variant.gt_types[index] 

                AF = AF/two_n
                AF_tuskless = AF_tuskless/two_n_tuskless
                AF_tusked = AF_tusked/two_n_tusked
                # is site polymorphic is tuskless, rare in tusked
                if AF_tuskless >= minAF_tuskless and AF_tuskless <= maxAF_tuskless:
                    if AF_tusked >= minAF_tusked and AF_tusked <= maxAF_tusked:
                        scaff = variant.CHROM
                        pos = variant.end
                        Depths = []
                        Genos = []
                        for index in range(len(variant.gt_types)):
                            if index2Group[index] == 0:
                                Genos.append(variant.gt_types[index])
                                Depths.append(variant.gt_depths[index])

                        positions[scaff].append(pos)
                        depthsDict[scaff].append(mean(Depths))

                        # feed  genotypes to functions to calculate HWE by phenotype
                        dev, het, hom = CalcHWEdev(Genos)
                        HWEdev[scaff].append(dev)
                        HETdev[scaff].append(het)
                        HOMdev[scaff].append(hom) 

    # assuming dict[scaffold] = [ list of values for scaffold ]
    # creates newDict[scaffold][window] = mean(values in winSize)
    HWEdevWin = makeWindowedDict(HWEdev, winSize, offset)
    HETdevWin = makeWindowedDict(HETdev, winSize, offset)
    HOMdevWin = makeWindowedDict(HOMdev, winSize, offset)
    depthsWin = makeWindowedDict(depthsDict, winSize, offset)

    positionsWin = defaultdict(lambda: defaultdict(tuple)) 
    for scaff in positions:
        # Dict[scaff] should be a (long) list that we will now dice up into windows
        wins = window(positions[scaff], winSize, offset) 
        cnt = 1 # label windows numerically
        for win in wins:
            positionsWin[scaff][cnt] = (win[0], win[-1])
            cnt += 1
        if scaff not in positionsWin:
            positionsWin[scaff][1] = ("NA", "NA")

    o = open(f"HWEdeviations_TusklessPrivateAlleles_{winSize}SNPwin_{offset}stepSize_minAF{minAF_tuskless}_maxAF{maxAF_tuskless}.txt", 'w')
    print("scaffold", end="\t", file=o) #1
    print("window", end="\t", file=o) #2
    print("lowPos", end="\t", file=o) #3
    print("highPos", end="\t", file=o) #4
    print("HWEdev", end="\t", file=o) #6
    print("HETdev", end="\t", file=o) #7
    print("HOMdev", end="\t", file=o) #8
    print("Depths", end="\n", file=o) #10

    for scaff in sorted(HWEdevWin):
        for win in sorted(HWEdevWin[scaff]):
            lower = positionsWin[scaff][win][0]
            higher = positionsWin[scaff][win][1]
            print(scaff, win, lower, higher, sep="\t", end="\t", file=o)
            print(HWEdevWin[scaff][win], HETdevWin[scaff][win], HOMdevWin[scaff][win], sep="\t", end="\t", file=o)
            print(depthsWin[scaff][win], sep="\t", file=o)
    o.close()

if __name__ == '__main__':
  main()
