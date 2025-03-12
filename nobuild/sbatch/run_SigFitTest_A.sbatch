#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --array=1-999
#SBATCH --cpus-per-task=8
#SBATCH --job-name=sigfittest
#SBATCH --mem=20GB
#SBATCH --partition EPYC
#SBATCH --output=/orfeo/cephfs/scratch/cdslab/ebusca00/GitHub/simbascule/nobuild/sbatch/logs/SigFitTest_A-%A-%a.log

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'

module load R/4.4.1

export PATH=~/miniconda3/bin:$PATH
conda activate basilica-env

Rscript /orfeo/cephfs/scratch/cdslab/ebusca00/GitHub/simbascule/nobuild/run_scripts/run_simuls_SigFitTest.R \
  $SLURM_ARRAY_TASK_ID 0
