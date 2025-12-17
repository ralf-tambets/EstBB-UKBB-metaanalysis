#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="MVP Coloc"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --array=0-13

module load python/3.12.3
source /path/to/venv/bin/activate

mkdir MVP_coloc

commands=(
"gpu-coloc -r --dir1 MVP/EUR_parquets --dir2 EST_UK_META/EstBB_parquets --results MVP_coloc/MVP_EUR_EstBB_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/EUR_parquets --dir2 EST_UK_META/UKBB_EUR_parquets --results MVP_coloc/MVP_EUR_UKBB_EUR_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/AFR_parquets --dir2 EST_UK_META/UKBB_AFR_parquets --results MVP_coloc/MVP_AFR_UKBB_AFR_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/AMR_parquets --dir2 EST_UK_META/UKBB_AMR_parquets --results MVP_coloc/MVP_AMR_UKBB_AMR_results.tsv --p12 1e-6 --H4 0"

"gpu-coloc -r --dir1 MVP/EUR_parquets --dir2 EST_UK_META/meta_EUR_parquets --results MVP_coloc/MVP_EUR_meta_EUR_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/AFR_parquets --dir2 EST_UK_META/meta_EUR_parquets --results MVP_coloc/MVP_AFR_meta_EUR_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/AMR_parquets --dir2 EST_UK_META/meta_EUR_parquets --results MVP_coloc/MVP_AMR_meta_EUR_results.tsv --p12 1e-6 --H4 0"

"gpu-coloc -r --dir1 MVP/META_parquets --dir2 EST_UK_META/EstBB_parquets --results MVP_coloc/MVP_META_EstBB_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/META_parquets --dir2 EST_UK_META/UKBB_EAS_parquets --results MVP_coloc/MVP_META_UKBB_EAS_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/META_parquets --dir2 EST_UK_META/UKBB_CSA_parquets --results MVP_coloc/MVP_META_UKBB_CSA_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/META_parquets --dir2 EST_UK_META/UKBB_EUR_parquets --results MVP_coloc/MVP_META_UKBB_EUR_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/META_parquets --dir2 EST_UK_META/UKBB_AFR_parquets --results MVP_coloc/MVP_META_UKBB_AFR_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/META_parquets --dir2 EST_UK_META/UKBB_AMR_parquets --results MVP_coloc/MVP_META_UKBB_AMR_results.tsv --p12 1e-6 --H4 0"
"gpu-coloc -r --dir1 MVP/META_parquets --dir2 EST_UK_META/meta_EUR_parquets --results MVP_coloc/MVP_META_meta_EUR_results.tsv --p12 1e-6 --H4 0"
)

eval "${commands[$SLURM_ARRAY_TASK_ID]}"
