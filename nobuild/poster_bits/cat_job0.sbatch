#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=EPYC
#SBATCH --job-name=0cat_array
#SBATCH --mem=50GB
#SBATCH --output=1cat_array-%A-%a.log

module load architecture/x86
module load conda
module load R/4.2.2

Rscript /u/cdslab/ebusca00/scratch_shared/basilica_pkgs/simbasilica/nobuild/poster_bits/run_newmodel_1606_wholeCat.R 0
