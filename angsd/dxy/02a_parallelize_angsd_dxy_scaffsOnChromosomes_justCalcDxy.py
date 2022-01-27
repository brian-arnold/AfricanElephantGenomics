#!/usr/bin/python -tt

import re
import sys
import os
sys.path.insert(1, "/tigress/bjarnold/AfricanElephants/scripts")
import functions as shFn
import time


def main():

    outDir = "./"
    bamListAll = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.All"
    bamListTuskless = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.MZtuskless"
    bamListTusked = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.MZtusked"
    ref = "/tigress/bjarnold/AfricanElephants/genome/Loxodonta_africana.loxAfr3.dna.toplevel.fa"
    angsd_dir = "/home/bjarnold/programs/angsd/angsd"
    script_dir = "/tigress/bjarnold/AfricanElephants/scripts/angsd/dxy/slidingWin"

    tasks = 1  # total number of tasks across all nodes  
    cpuPerTask = 1 # cpu-cores per task (>1 if multi-threaded tasks)  
    t = "05:00:00" 
    mem = 10000

    scaffs = []
    ih = open("/tigress/bjarnold/AfricanElephants/genome/loxAfr4/ScaffoldsInChromosomes.txt",'r')
    for line in ih:
        if not line.startswith("scaffoldName"):
            line = line.split()
            scaffs.append(line[0])
    ih.close()

    #cnt = 1
    for scaffold in scaffs:
        #if cnt == 1:
        #f = "Dxy_" + scaffold + "_winSize10000.txt"
        #if not os.path.isfile(f):
        command = ""
        command = command + "python3 " + script_dir + "/02_dxy_slidingWindow_SNPwin.py "
        command = command + scaffold + ".pos.gz " + scaffold + "_tuskless.mafs.gz " + scaffold + "_tusked.mafs.gz\n"
        print(scaffold)
        shFn.create_job_script(scaffold, outDir, tasks, cpuPerTask, t, mem, command) 
        time.sleep(3)
        #cnt += 1

if __name__ == '__main__':
  main()
