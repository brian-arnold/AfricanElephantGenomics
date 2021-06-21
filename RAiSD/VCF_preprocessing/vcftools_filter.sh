#!/bin/bash
#SBATCH -J template
#SBATCH -o out
#SBATCH -e err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time 1-00:00:00        # DAYS-HOURS:MINUTES:SECONDS


source /home/bjarnold/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

vcftools \
--vcf /scratch/gpfs/bjarnold/AE/VCF/Final_refiltered_annotated.vcf \
--remove-filtered-all \
--minDP 8 \
--maxDP 80 \
--indv 0045B \
--indv 2983B \
--indv 2984B \
--indv G17A \
--indv G19A \
--indv G22A \
--indv T2B \
--out tusklessMZ_DP8_2 \
--recode
# NOTE: exluded 0045B b/c sib of 2986A which had less miss data

#--minGQ 15 \


