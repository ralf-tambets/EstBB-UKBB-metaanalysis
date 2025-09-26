#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="FinnGen coloc"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --array=0-6

datasets=(EstBB UKBB_AFR UKBB_AMR UKBB_CSA UKBB_EAS UKBB_EUR UKBB_MID)

current_dataset=${datasets[$SLURM_ARRAY_TASK_ID]}

module load python/3.12.3
source /path/to/venv/bin/activate # requires pip install gpu-coloc
mkdir FinnGen_lbf_coloc
gpu-coloc -r --dir1 /path/to/FinnGen_lbf_parquets --dir2 EST_UK_META/${current_dataset}_parquets --results FinnGen_lbf_coloc/${current_dataset}_FinnGen_lbfs_results.tsv --p12 1e-6 --H4 0
mkdir FinnGen_abf_coloc
gpu-coloc -r --dir1 /path/to/FinnGen_abf_parquets --dir2 EST_UK_META/${current_dataset}_parquets --results FinnGen_abf_coloc/${current_dataset}_FinnGen_abfs_results.tsv --p12 1e-6 --H4 0

