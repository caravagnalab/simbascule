#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --array=0-8
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:8
#SBATCH --job-name=1dn_array
#SBATCH --mem=50GB
#SBATCH --mail-type=ALL
#SBATCH --output=1dn_array-%A-%a.log

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
module load R/4.2.2

# Rscript /u/cdslab/ebusca00/scratch_shared/basilica_pkgs/simbasilica/nobuild/script_test/run_newmodel_3105.R
