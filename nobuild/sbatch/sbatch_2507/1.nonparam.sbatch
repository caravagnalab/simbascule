#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --array=0-11
#SBATCH --partition=cpuq
#SBATCH --job-name=1.nonparam
#SBATCH --mem=20GB
#SBATCH --output=/home/elena.buscaroli/GitHub/simbasilica/nobuild/sbatch_2507/logs_ht/1.nonparam-%A-%a.log

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'

# module load cuda11.1/toolkit/11.1.1
module load R/4.2.3

export PATH=~/miniconda3/bin:$PATH
conda activate basilica-env

Rscript /home/elena.buscaroli/GitHub/simbasilica/nobuild/run_scripts/run_power_test.R \
  $SLURM_ARRAY_TASK_ID \
  nonparametric \
  flat_clust \
  2507

