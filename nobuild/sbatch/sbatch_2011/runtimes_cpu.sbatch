#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1
#SBATCH --array=1-9
#SBATCH --mem=50G
#SBATCH --job-name=runtimes
#SBATCH --output=~/GitHub/simbasilica/nobuild/sbatch/logs/runtimes-%A-%a.log

module load R/4.2.3

Rscript ~/GitHub/simbasilica/nobuild/run_scripts/run_execution_time.R \
  $SLURM_ARRAY_TASK_ID \
  runtimes

