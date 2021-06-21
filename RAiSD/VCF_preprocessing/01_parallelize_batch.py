#!/usr/bin/python -tt

import re
import sys
import os
import time
sys.path.insert(1, "/tigress/bjarnold/AfricanElephants/scripts")
import functions as shFn

def main():

    tasks = 1  # total number of tasks across all nodes
    cpuPerTask = 1 # cpu-cores per task (>1 if multi-threaded tasks)
    t = "08:00:00"
    mem = 30000

    outDir = "."
    inFile = "/scratch/gpfs/bjarnold/AE/workflows/RAiSD/tusklessMZ_DP8.recode.vcf"
    scaffolds = []
    f = open("/tigress/bjarnold/AfricanElephants/genome/scaffolds.txt",'r')
    for line in f:
        line = line.split()
        scaff = line[0]
        scaff = scaff.strip(">")
        scaffolds.append(scaff)
    count = 1
    command = "source /home/bjarnold/miniconda3/etc/profile.d/conda.sh\n"
    command = command + "conda activate vcftools\n"
    for scaff in scaffolds:
        if count >= 10:
            command = command + "vcftools --vcf " + inFile + " "
            command = command + "--chr " + scaff + " "
            command = command + "--out " + scaff + " "
            command = command + "--recode\n"
            shFn.create_job_script(scaff, outDir, tasks, cpuPerTask, t, mem, command)	
            # re-initialize command
            count = 1
            command = "source /home/bjarnold/miniconda3/etc/profile.d/conda.sh\n"
            command = command + "conda activate vcftools\n"
        else:
            command = command + "vcftools --vcf " + inFile + " "
            command = command + "--chr " + scaff + " "
            command = command + "--out " + scaff + " "
            command = command + "--recode\n"
            count += 1
    # submit straglers
    shFn.create_job_script(scaff, outDir, tasks, cpuPerTask, t, mem, command)	


if __name__ == '__main__':
  main()
