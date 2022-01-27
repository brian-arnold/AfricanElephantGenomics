#!/usr/bin/python -tt

import re
import sys
import time
import os
sys.path.insert(1, "/tigress/bjarnold/AfricanElephants/scripts")
import functions as shFn


def main():

    outDir = "./"
    bamListAll = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.All"
    bamListTuskless = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.MZtuskless"
    bamListTusked = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.MZtusked"
    ref = "/tigress/bjarnold/AfricanElephants/genome/Loxodonta_africana.loxAfr3.dna.toplevel.fa"
    angsd_dir = "/home/bjarnold/programs/angsd/angsd"
    script_dir = "/tigress/bjarnold/AfricanElephants/scripts/angsd/dxy"

    tasks = 1  # total number of tasks across all nodes  
    cpuPerTask = 2 # cpu-cores per task (>1 if multi-threaded tasks)  
    t = "08:00:00" 
    mem = 30000

    scaffsChrom = []
    ih = open("/tigress/bjarnold/AfricanElephants/genome/loxAfr4/ScaffoldsInChromosomes.txt",'r')
    for line in ih:
        if not line.startswith("scaffoldName"):
            line = line.split()
            scaffsChrom.append(line[0])
    ih.close()

    scaffsTotal = []
    ih = open("/tigress/bjarnold/AfricanElephants/genome/scaffolds.txt",'r')
    for line in ih:
            line = line.split()
            s = line[0].strip(">")
            scaffsTotal.append(s)
    ih.close()

    scaffsOrphan = []
    for scaff in scaffsTotal:
        if scaff not in scaffsChrom:
            scaffsOrphan.append(scaff)

    for scaffold in scaffsOrphan:
        #f = "Dxy_" + scaffold + "_winSize10000.txt"
        #if not os.path.isfile(f):
        #if os.path.getsize(f) == 0:
        command = ""
        command = command + "ulimit -s $((2**20))\n"
        command = command + "set -o errexit\n"
        command = command + "ulimit -c unlimited\n"

        # get allele frequencies across entire population
        command = command + f"{angsd_dir}/angsd -bam {bamListAll} "
        command = command + "-doMajorMinor 4 "
        command = command + "-doMaf 2 "
        command = command + "-GL 1 "
        command = command + f"-ref {ref} "
        command = command + "-doCounts 1 "
        command = command + "-minMaf 0.1 "
        command = command + "-setMinDepth 56 "
        #command = command + "-minInd 13 "
        #command = command + "-minIndDepth 2 "
        command = command + "-minQ 20 "
        command = command + "-minMapQ 20 "
        command = command + "-remove_bads "
        command = command + "-SNP_pval 1e-6 "
        command = command + "-r " + scaffold + " "
        command = command + "-P " + str(cpuPerTask) + " "
        command = command + "-out " + scaffold + "_All\n"

        # index sites for using with single populations
        command = command + "python3 " + script_dir + "/01b_reformat_maf_2_sites.py " + scaffold + "\n"
        command = command + "sleep 30\n"
        command = command + f"{angsd_dir}/angsd sites index " + scaffold + ".sites\n"

        # get callable sites
        command = command + f"{angsd_dir}/angsd -bam " + bamListAll + " "
        command = command + "-doCounts 1 "
        command = command + "-dumpCounts 1 "
        #command = command + "-minMaf 0.1 "
        #command = command + "-minInd 13 "
        #command = command + "-minIndDepth 2 "
        command = command + "-minQ 20 "
        command = command + "-minMapQ 20 "
        command = command + "-remove_bads "
        command = command + "-r " + scaffold + " "
        command = command + "-P " + str(cpuPerTask) + " "
        command = command + "-out " + scaffold + "\n"

        # get callable sites Tuskless
        command = command + f"{angsd_dir}/angsd -bam " + bamListTuskless + " "
        command = command + "-doCounts 1 "
        command = command + "-dumpCounts 1 "
        #command = command + "-minMaf 0.1 "
        #command = command + "-minInd 13 "
        #command = command + "-minIndDepth 2 "
        command = command + "-minQ 20 "
        command = command + "-minMapQ 20 "
        command = command + "-remove_bads "
        command = command + "-r " + scaffold + " "
        command = command + "-P " + str(cpuPerTask) + " "
        command = command + "-out " + scaffold + "_tuskless\n"

        # get callable sites Tusked
        command = command + f"{angsd_dir}/angsd -bam " + bamListTusked + " "
        command = command + "-doCounts 1 "
        command = command + "-dumpCounts 1 "
        #command = command + "-minMaf 0.1 "
        #command = command + "-minInd 13 "
        #command = command + "-minIndDepth 2 "
        command = command + "-minQ 20 "
        command = command + "-minMapQ 20 "
        command = command + "-remove_bads "
        command = command + "-r " + scaffold + " "
        command = command + "-P " + str(cpuPerTask) + " "
        command = command + "-out " + scaffold + "_tusked\n"

        # call sites for tuskless
        command = command + f"{angsd_dir}/angsd -bam " + bamListTuskless + " "
        command = command + "-doMajorMinor 4 "
        command = command + f"-ref {ref} "
        command = command + "-doMaf 2 "
        command = command + "-GL 1 "
        command = command + "-r " + scaffold + " "
        command = command + "-sites " + scaffold + ".sites "
        command = command + "-P " + str(cpuPerTask) + " "
        command = command + "-out " + scaffold + "_tuskless\n"

        # call sites for tusked 
        command = command + f"{angsd_dir}/angsd -bam " + bamListTusked + " "
        command = command + "-doMajorMinor 4 "
        command = command + f"-ref {ref} "
        command = command + "-doMaf 2 "
        command = command + "-GL 1 "
        command = command + "-r " + scaffold + " "
        command = command + "-sites " + scaffold + ".sites "
        command = command + "-P " + str(cpuPerTask) + " "
        command = command + "-out " + scaffold + "_tusked\n"

        #command = command + "python3 " + script_dir + "/02_dxy.py "
        command = command + scaffold + ".pos.gz " + scaffold + "_tuskless.mafs.gz " + scaffold + "_tusked.mafs.gz\n" 

        print(scaffold)
        shFn.create_job_script(scaffold, outDir, tasks, cpuPerTask, t, mem, command) 
        #time.sleep(20)

if __name__ == '__main__':
  main()
