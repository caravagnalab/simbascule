#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --array=4-6
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:3
#SBATCH --job-name=2dn_array
#SBATCH --mem=50GB
#SBATCH --output=/home/elena.buscaroli/GitHub/simbasilica/nobuild/poster_bits/logs_ht/2dn_array-%A-%a.log

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'

module load cuda11.1/toolkit/11.1.1
module load R/4.2.3

export PATH=~/miniconda3/bin:$PATH
conda activate basilica-env

Rscript /home/elena.buscaroli/GitHub/simbasilica/nobuild/poster_bits/run_newmodel_2606_dn.R  $SLURM_ARRAY_TASK_ID
