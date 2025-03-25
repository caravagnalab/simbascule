#!/bin/bash
#SBATCH --job-name=mutsigR
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=10
#SBATCH --time=60:00:00
#SBATCH --output=signaturetools_real.out
#SBATCH --mem=30G

cd $SLURM_SUBMIT_DIR

source /u/area/evillegas/.bashrc
#module load conda
conda init
conda activate rmutsig

echo $(date)
Rscript cosmic_fit_real.R > logs/cosmic_output_real.txt 2>logs/cosmic_error_real.txt
# Rscript fitms_real.R > logs/output_fitms.txt 2>logs/error_fitms.txt
echo $(date)
