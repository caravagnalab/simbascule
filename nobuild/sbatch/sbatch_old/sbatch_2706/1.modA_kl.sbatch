#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --array=0-2
#SBATCH --partition=cpuq
#SBATCH --job-name=1.A_cos
#SBATCH --mem=50GB
#SBATCH --output=/home/elena.buscaroli/GitHub/simbascule/nobuild/sbatch_2706/logs_ht/1.A_cos-%A-%a.log

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
conda activate bascule-env

Rscript /home/elena.buscaroli/GitHub/simbascule/nobuild/run_scripts/run_dn.R $SLURM_ARRAY_TASK_ID new_hier KL hier_clust


