#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="FinnGen+MVP+UKBB meta Coloc"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

module load python/3.12.3
source /path/to/venv/bin/activate

mkdir FinnGen+MVP+UKBB_coloc
gpu-coloc -r --dir1 FinnGen+MVP+UKBB/FinnGen+MVP+UKBB_parquets --dir2 EST_UK_META/meta_EUR_parquets --results FinnGen+MVP+UKBB_coloc/meta_EUR_parquets_FinnGen+MVP+UKBB_META_results.tsv --p12 1e-6 --H4 0
