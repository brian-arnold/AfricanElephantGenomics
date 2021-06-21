#!/usr/bin/python -tt

import re
import sys
import os
import time
sys.path.insert(1, "/tigress/bjarnold/AfricanElephants/scripts")
import functions as shFn


def main():

    outDir = "./"
    bamListTuskless = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.MZtuskless"
    bamListTusked = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.MZtusked"
    ref = "/tigress/bjarnold/AfricanElephants/genome/Loxodonta_africana.loxAfr3.dna.toplevel.fa"
    angsd_dir = "/home/bjarnold/programs/angsd/angsd"

    tasks = 1  # total number of tasks across all nodes  
    cpuPerTask = 1 # cpu-cores per task (>1 if multi-threaded tasks)  
    t = "05:00:00" 
    #mem = 10000
    mem = 40000

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
        f = scaffold + ".fst10kbwin2kbslide.txt"
        #if scaffold != "scaffold_17":
        if not os.path.isfile(f):
            # exanple taken from http://www.popgen.dk/angsd/index.php/Fst
            command = ""
            command = command + f"{angsd_dir}/angsd -bam {bamListTuskless} "
            command = command + f"-anc {ref} "
            command = command + f"-doSaf 1 "
            command = command + f"-GL 1 "
            command = command + f"-r {scaffold} "
            command = command + f"-out {scaffold}_tuskless\n"

            command = command + f"{angsd_dir}/angsd -bam {bamListTusked} "
            command = command + f"-anc {ref} "
            command = command + f"-doSaf 1 "
            command = command + f"-GL 1 "
            command = command + f"-r {scaffold} "
            command = command + f"-out {scaffold}_tusked\n"

            command = command + f"{angsd_dir}/misc/realSFS " 
            command = command + f"{scaffold}_tuskless.saf.idx "
            command = command + f"{scaffold}_tusked.saf.idx "
            command = command + f"> {scaffold}tuskless.tusked.ml\n"

            command = command + f"{angsd_dir}/misc/realSFS fst index "
            command = command + f"{scaffold}_tuskless.saf.idx " 
            command = command + f"{scaffold}_tusked.saf.idx " 
            command = command + f"-r {scaffold} " 
            command = command + f"-sfs {scaffold}tuskless.tusked.ml " 
            command = command + f"-fstout {scaffold}\n"

            command = command + f"{angsd_dir}/misc/realSFS fst stats2 "
            command = command + f"{scaffold}.fst.idx " 
            command = command + f"-win 10000 " 
            command = command + f"-step 2000 " 
            command = command + f"> {scaffold}.fst10kbwin2kbslide.txt\n" 

            print(scaffold)
            shFn.create_job_script(scaffold, outDir, tasks, cpuPerTask, t, mem, command) 
            #time.sleep(5)

if __name__ == '__main__':
  main()
