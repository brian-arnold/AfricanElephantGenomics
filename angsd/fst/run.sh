#!/bin/bash
#SBATCH -J angsd
#SBATCH -o out
#SBATCH -e err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time 1-00:00:00        # DAYS-HOURS:MINUTES:SECONDS

# run this in the directory where you want output files
python3 /tigress/bjarnold/AfricanElephants/scripts/angsd/fst/01a_parallelize_angsd_fst_scaffsOnChromosomes.py
python3 /tigress/bjarnold/AfricanElephants/scripts/angsd/fst/01b_parallelize_angsd_fst_scaffsOrphans.py


