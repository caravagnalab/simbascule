#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=2
#SBATCH --cpus-per-task=1
#SBATCH --array=0-4
#SBATCH --partition=CPU
#SBATCH --job-name=1cat_array
#SBATCH --mem=50GB
#SBATCH --output=1cat_array-%A-%a.log

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'

module load architecture/Intel
module load cuda
module load R/4.2.2

Rscript /u/cdslab/ebusca00/scratch_shared/basilica_pkgs/simbasilica/nobuild/poster_bits/run_newmodel_1606_wholeCat.R  $SLURM_ARRAY_TASK_ID
