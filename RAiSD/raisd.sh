#!/bin/bash
#SBATCH -J raisd
#SBATCH -o out
#SBATCH -e err
#SBATCH -p holy-info
#SBATCH -n 1
#SBATCH -t 5000 
#SBATCH --mem=5000

python3 01_parallelize.py

