#!/bin/bash
#SBATCH --time=30:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:1
#SBATCH --job-name=dmm_real
#SBATCH --mem=30GB
#SBATCH --output=/home/elena.buscaroli/GitHub/simbascule/nobuild/sbatch/sbatch_1809/logs_ht/dmm_real-%A.log

module load cuda11.1/toolkit/11.1.1
module load R/4.2.3

export PATH=~/miniconda3/bin:$PATH
conda activate bascule-env

Rscript /home/elena.buscaroli/GitHub/simbascule/nobuild/run_scripts/run_real_data.sbs.R \
  -1 \
  clust \
  dmm.2809
