#!/usr/bin/python -tt

import re
import sys
import os
import functions as shFn

def main():

	queue = "holy-info,shared"
	#queue = "bos-info,shared"
	time = 80
	mem = 5000
	nCores = 1
	outDir = "."
	#inFile = "/n/holyscratch01/informatics/bjarnold/AfricanElephants/workflow/VCFs/4thRun_AddedMZind/VCFs_for_RAiSD/AllMZ.recode.vcf"
	#inFile = "/n/holyscratch01/informatics/bjarnold/AfricanElephants/workflow/VCFs/4thRun_AddedMZind/VCFs_for_RAiSD/tuskedMZ.recode.vcf"
	#inFile = "/n/holylfs/LABS/informatics/bjarnold/AfricanElephants/workflow/VCFs/4thRun_AddedMZind/VCFs_for_RAiSD/tusklessMZ.recode.vcf"
	inFile = "/n/holylfs/LABS/informatics/bjarnold/AfricanElephants/workflow/VCFs/4thRun_AddedMZind/VCFs_for_RAiSD/tusklessMZ_8X.recode.vcf"
	scaffolds = []
	f = open("/n/holylfs/LABS/informatics/bjarnold/AfricanElephants/genome/scaffolds.txt",'r')
	for line in f:
            line = line.split()
            scaff = line[0]
            scaff = scaff.strip(">")
            scaffolds.append(scaff)
	count = 0
	for scaff in scaffolds:
            count +=1
            #if count == 1:
            outFile = "./tusklessMZ_splitByScaff_8X/" + scaff + ".recode.vcf"
            if not os.path.isfile(outFile):
                command = "module load vcftools/0.1.14-fasrc01\n"
                command = command + "vcftools --vcf " + inFile + " "
                command = command + "--chr " + scaff + " "
                command = command + "--out " + scaff + " "
                command = command + "--recode\n"

                print(outFile)
                print(scaff)
                #shFn.create_job_script(scaff, outDir, queue, nCores, time, mem, command)	

if __name__ == '__main__':
  main()
