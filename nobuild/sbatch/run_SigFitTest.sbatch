#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --array=1-1620
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sigfittest
#SBATCH --mem=20GB
#SBATCH --partition boost_usr_prod
#SBATCH --output=/leonardo_scratch/fast/uTS25_Caravagn/ebuscaro/GitHub/simbascule/nobuild/sbatch/logs/SigFitTest-%A-%a.log

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'


module load anaconda3

source activate bascule-env

module load spack
module load r


Rscript /leonardo_scratch/fast/uTS25_Caravagn/ebuscaro/GitHub/simbascule/nobuild/run_scripts/run_simuls_SigFitTest.R \
  $SLURM_ARRAY_TASK_ID 0

