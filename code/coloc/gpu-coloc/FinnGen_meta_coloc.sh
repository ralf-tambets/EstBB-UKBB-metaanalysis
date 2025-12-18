#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="FinnGen meta_EUR coloc"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load python/3.12.3
source /path/to/venv/bin/activate # requires pip install gpu-coloc
mkdir FinnGen_lbf_coloc
gpu-coloc -r --dir1 /path/to/FinnGen_lbf_parquets --dir2 EST_UK_META/meta_EUR_parquets --results FinnGen_lbf_coloc/meta_EUR_FinnGen_lbfs_results.tsv --p12 1e-6 --H4 0
mkdir FinnGen_abf_coloc
gpu-coloc -r --dir1 /path/to/FinnGen_abf_parquets --dir2 EST_UK_META/meta_EUR_parquets --results FinnGen_abf_coloc/meta_EUR_FinnGen_abfs_results.tsv --p12 1e-6 --H4 0

