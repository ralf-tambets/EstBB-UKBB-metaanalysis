#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="UKBPPP coloc"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --array=0-20

mkdir UKBPPP_coloc

module load python/3.12.3
source /path/to/venv/bin/activate

commands=(
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/European_(discovery)_parquets\" --dir2 EST_UK_META/UKBB_EUR_parquets --results UKBPPP_coloc/UKBPPP_EUR_UKBB_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/European_(discovery)_parquets\" --dir2 EST_UK_META/EstBB_parquets --results UKBPPP_coloc/UKBPPP_EUR_EstBB_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/African_parquets --dir2 EST_UK_META/UKBB_AFR_parquets --results test_UKBPPP_AFR_UKBB_AFR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/Middle_East_parquets --dir2 EST_UK_META/UKBB_MID_parquets --results UKBPPP_coloc/UKBPPP_MID_UKBB_MID_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/American_parquets --dir2 EST_UK_META/UKBB_AMR_parquets --results UKBPPP_coloc/UKBPPP_AMR_UKBB_AMR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/Central_South_Asian_parquets --dir2 EST_UK_META/UKBB_CSA_parquets --results UKBPPP_coloc/UKBPPP_CSA_UKBB_CSA_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/Eeas_Asian_parquets --dir2 EST_UK_META/UKBB_EAS_parquets --results UKBPPP_coloc/UKBPPP_EAS_UKBB_EAS_results.tsv --p12 1e-6 --H4 0.8"

  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/Combined_parquets\" --dir2 EST_UK_META/UKBB_EUR_parquets --results UKBPPP_coloc/UKBPPP_Combined_UKBB_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/Combined_parquets\" --dir2 EST_UK_META/EstBB_parquets --results UKBPPP_coloc/UKBPPP_Combined_EstBB_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/Combined_parquets\" --dir2 EST_UK_META/UKBB_AFR_parquets --results UKBPPP_coloc/UKBPPP_Combined_UKBB_AFR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/Combined_parquets\" --dir2 EST_UK_META/UKBB_MID_parquets --results UKBPPP_coloc/UKBPPP_Combined_UKBB_MID_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/Combined_parquets\" --dir2 EST_UK_META/UKBB_AMR_parquets --results UKBPPP_coloc/UKBPPP_Combined_UKBB_AMR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/Combined_parquets\" --dir2 EST_UK_META/UKBB_CSA_parquets --results UKBPPP_coloc/UKBPPP_Combined_UKBB_CSA_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/Combined_parquets\" --dir2 EST_UK_META/UKBB_EAS_parquets --results UKBPPP_coloc/UKBPPP_Combined_UKBB_EAS_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/Combined_parquets\" --dir2 EST_UK_META/meta_EUR_parquets --results UKBPPP_coloc/UKBPPP_Combined_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"

  "gpu-coloc -r --dir1 \"/path/to/ukbppp_parquets_fixed/European_(discovery)_parquets\" --dir2 EST_UK_META/meta_EUR_parquets --results UKBPPP_coloc/UKBPPP_EUR_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/African_parquets --dir2 EST_UK_META/meta_EUR_parquets --results UKBPPP_coloc/UKBPPP_AFR_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/Middle_East_parquets --dir2 EST_UK_META/meta_EUR_parquets --results UKBPPP_coloc/UKBPPP_MID_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/American_parquets --dir2 EST_UK_META/meta_EUR_parquets --results UKBPPP_coloc/UKBPPP_AMR_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/Central_South_Asian_parquets --dir2 EST_UK_META/meta_EUR_parquets --results UKBPPP_coloc/UKBPPP_CSA_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/ukbppp_parquets_fixed/Eeas_Asian_parquets --dir2 EST_UK_META/meta_EUR_parquets --results UKBPPP_coloc/UKBPPP_EAS_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
)

eval "${commands[$SLURM_ARRAY_TASK_ID]}"
