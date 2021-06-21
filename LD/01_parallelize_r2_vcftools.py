#!/usr/bin/python -tt

import re
import sys
import os
sys.path.insert(1, "/tigress/bjarnold/AfricanElephants/scripts")
import functions as shFn


def main():

    outDir = "./"
    tasks = 1  # total number of tasks across all nodes  
    cpuPerTask = 1 # cpu-cores per task (>1 if multi-threaded tasks)  
    time = "12:00:00" 
    mem = 40000

    vcf = "/scratch/gpfs/bjarnold/AE/VCF/Final_refiltered_annotated.vcf"

    tuskless = ["0045B", "2983B", "2984B", "G17A", "G19A", "G22A", "T2B"]
    tusked = ["2981B", "2982B", "2985B", "G18A", "G20A", "G21A"]
    samples = tuskless + tusked

    scaffs = []
    ih = open("/tigress/bjarnold/AfricanElephants/genome/loxAfr4/ScaffoldsInChromosomes.txt",'r')
    for line in ih:
        if not line.startswith("scaffoldName"):
            line = line.split()
            scaffs.append(line[0])
    ih.close()
    count = 0
    for scaffold in scaffs:
        count += 1
        if count >= 1:
            command = ""
            command = command + "source /home/bjarnold/miniconda3/etc/profile.d/conda.sh\n"
            command = command + "conda activate vcftools\n"
            command = command + "vcftools "
            command = command + f"--vcf {vcf} "
            command = command + "--remove-indels "
            for i in samples:
                command = command + f"--indv {i} "
            command = command + "--minDP 8 "
            command = command + "--maxDP 80 "
            command = command + "--maf 0.05 "
            command = command + "--max-maf 0.95 "
            command = command + "--max-missing-count 0 "
            command = command + "--ld-window-bp 200000 "
            command = command + "--remove-filtered-all "
            command = command + f"--chr {scaffold} "
            command = command + f"--out {scaffold} "
            command = command + "--geno-r2 "

            shFn.create_job_script(scaffold, outDir, tasks, cpuPerTask, time, mem, command) 

if __name__ == '__main__':
  main()
