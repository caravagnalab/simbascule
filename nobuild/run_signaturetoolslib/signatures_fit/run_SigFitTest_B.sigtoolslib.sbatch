#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --array=0-620
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sigfittest
#SBATCH --mem=20GB
#SBATCH --partition Main,lovelace
#SBATCH --output=/share/cdslab/elenab/GitHub/simbascule/nobuild/run_signaturetoolslib/logs/SigFitTest_B-%A-%a.log

## Leonardo params
# --partition boost_usr_prod
# --output=/leonardo_scratch/fast/uTS25_Caravagn/ebuscaro/GitHub/simbascule/nobuild/sbatch/logs.FAST/SigFitTest_B-%A-%a.log
# serial partition: lrd_all_serial; --gres=tmpfs:1g

echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'

## Leonardo modules
# module load anaconda3/
# source activate bascule-env
# module load spack
# module load r

## Demetra modules
module load r/4.4.0-gcc-13.2.0-hymdfgk
conda activate bascule-env

## Leonardo script
# Rscript /leonardo_scratch/fast/uTS25_Caravagn/ebuscaro/GitHub/simbascule/nobuild/run_scripts/run_simuls_SigFitTest.FAST.R \

## Demetra script
Rscript /share/cdslab/elenab/GitHub/simbascule/nobuild/run_signaturetoolslib/run_sigtoolslib.R \
  $SLURM_ARRAY_TASK_ID 1000 SigFitTest

