#!/usr/bin/python -tt

import re
import sys
import time
import os
sys.path.insert(1, "/tigress/bjarnold/AfricanElephants/scripts")
import functions as shFn

def main():

    winSize = 50
    slack = 1

    tasks = 1  # total number of tasks across all nodes
    cpuPerTask = 1 # cpu-cores per task (>1 if multi-threaded tasks)
    t = "05:00:00"
    mem = 5000

    outDir = "."
    inDir = "/scratch/gpfs/bjarnold/AE/workflows/RAiSD/VCFs_splitByScaff/"
    scaffolds = []
    f = open("/tigress/bjarnold/AfricanElephants/genome/scaffolds.txt",'r')
    for line in f:
        line = line.split()
        scaff = line[0]
        scaff = scaff.strip(">")
        scaffolds.append(scaff)
    count = 0
    for scaff in scaffolds:
        #count +=1
        #if count == 1:
        command = "/home/bjarnold/programs/RAiSD/raisd-master/RAiSD "
        command = command + "-I " + inDir + scaff + ".recode.vcf "
        command = command + "-n " + scaff + " "
        command = command + "-w " + str(winSize) + " "
        command = command + "-c " + str(slack) + " "
        command = command + "-f -R -D\n"
        shFn.create_job_script(scaff, outDir, tasks, cpuPerTask, t, mem, command)	
        time.sleep(3)

if __name__ == '__main__':
  main()
