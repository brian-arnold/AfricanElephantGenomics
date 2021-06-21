#!/bin/bash
#SBATCH -J template
#SBATCH -o out
#SBATCH -e err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=30G         # memory per cpu-core (4G is default)
#SBATCH --time 1-00:00:00        # DAYS-HOURS:MINUTES:SECONDS

source /home/bjarnold/miniconda3/etc/profile.d/conda.sh
conda activate overlaps

python3 /tigress/bjarnold/AfricanElephants/scripts/overlap/01_combine_tests.py
python3 /tigress/bjarnold/AfricanElephants/scripts/overlap/02_overlap_tests.py 95
python3 /tigress/bjarnold/AfricanElephants/scripts/overlap/03_overlap_dental_add_chr.py



