#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --array=1-3
#SBATCH --cpus-per-task=8
#SBATCH --job-name=stl_GEL
#SBATCH --mem=20GB
#SBATCH --partition EPYC
#SBATCH --output=/orfeo/cephfs/scratch/cdslab/ebusca00/GitHub/simbascule/nobuild/run_signaturetoolslib/signatures_real/logs/GEL-%A-%a.log

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'


## Orfeo modules
module load R/4.4.1

## Orfeo script
Rscript ~/GitHub/simbascule/nobuild/run_signaturetoolslib/signatures_real/run_sigtoolslib_real.R \
  $SLURM_ARRAY_TASK_ID GEL

