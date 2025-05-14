#!/bin/bash
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --array=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sft_stl_TEST
#SBATCH --mem=10GB
#SBATCH --partition EPYC
#SBATCH --output=/orfeo/scratch/cdslab/ebusca00/GitHub/simbascule/nobuild/run_signaturetoolslib/signatures_extraction/logs/sft_stl_TEST-%A-%a.log

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

# ## Demetra modules
# module load r/4.4.0-gcc-13.2.0-hymdfgk
# conda activate bascule-env

## Orfeo modules
module load R/4.4.1


## Leonardo script
# Rscript /leonardo_scratch/fast/uTS25_Caravagn/ebuscaro/GitHub/simbascule/nobuild/run_scripts/run_simuls_SigFitTest.FAST.R \

# ## Demetra script
# Rscript /share/cdslab/elenab/GitHub/simbascule/nobuild/run_signaturetoolslib/run_sigtoolslib.R \
#   $SLURM_ARRAY_TASK_ID 0 SigFitTest

## Orfeo script
Rscript ~/GitHub/simbascule/nobuild/run_signaturetoolslib/signatures_extraction/run_sigtoolslib_extraction.R \
  $SLURM_ARRAY_TASK_ID 0 SigFitTest
