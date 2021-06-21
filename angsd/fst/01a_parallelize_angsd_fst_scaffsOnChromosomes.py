#!/usr/bin/python -tt

import re
import sys
import os
sys.path.insert(1, "/tigress/bjarnold/AfricanElephants/scripts")
import functions as shFn


def main():

    outDir = "./"
    bamListTuskless = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.MZtuskless"
    bamListTusked = "/tigress/bjarnold/AfricanElephants/scripts/angsd/bam.list.MZtusked"
    ref = "/tigress/bjarnold/AfricanElephants/genome/Loxodonta_africana.loxAfr3.dna.toplevel.fa"
    angsd_dir = "/home/bjarnold/programs/angsd/angsd"

    tasks = 1  # total number of tasks across all nodes  
    cpuPerTask = 10 # cpu-cores per task (>1 if multi-threaded tasks)  
    time = "23:00:00" 
    mem = 50000

    scaffs = []
    ih = open("/tigress/bjarnold/AfricanElephants/genome/loxAfr4/ScaffoldsInChromosomes.txt",'r')
    for line in ih:
        if not line.startswith("scaffoldName"):
            line = line.split()
            scaffs.append(line[0])
    ih.close()

    for scaffold in scaffs:
        f = scaffold + ".fst10kbwin2kbslide.txt"
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
            shFn.create_job_script(scaffold, outDir, tasks, cpuPerTask, time, mem, command) 

if __name__ == '__main__':
  main()
