#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=24
#SBATCH --partition=GPU
#SBATCH --job-name=2cat_job
#SBATCH --mem=50GB
#SBATCH --output=2cat_job-%A-%a.log

module load architecture/x86
module load conda
module load cuda
module load R/4.2.2

Rscript /u/cdslab/ebusca00/scratch_shared/basilica_pkgs/simbasilica/nobuild/poster_bits/run_newmodel_1606_wholeCat.R NULL

