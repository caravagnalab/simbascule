#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --array=1-9
#SBATCH --mem=50G
#SBATCH --job-name=runtimes_gpu
#SBATCH --output=/u/cdslab/ebusca00/GitHub/simbasilica/nobuild/sbatch/logs/runtimes_gpus-%A-%a.log

module load R/4.2.3
module load cuda

Rscript ~/GitHub/simbasilica/nobuild/run_scripts/run_execution_time.R \
  $SLURM_ARRAY_TASK_ID \
  runtimes \
  gpu

