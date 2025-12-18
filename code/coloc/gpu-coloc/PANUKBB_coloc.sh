#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="PANUKBB coloc"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --array=0-20

module load python/3.12.3
source path/to/venv/bin/activate

mkdir PANUKBB_coloc

commands=(
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/EUR_parquets\" --dir2 EST_UK_META/UKBB_EUR_parquets --results PANUKBB_coloc/PANUKBB_EUR_UKBB_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/EUR_parquets\" --dir2 EST_UK_META/EstBB_parquets --results PANUKBB_coloc/PANUKBB_EUR_EstBB_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/AFR_parquets --dir2 EST_UK_META/UKBB_AFR_parquets --results PANUKBB_AFR_UKBB_AFR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/MID_parquets --dir2 EST_UK_META/UKBB_MID_parquets --results PANUKBB_coloc/PANUKBB_MID_UKBB_MID_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/AMR_parquets --dir2 EST_UK_META/UKBB_AMR_parquets --results PANUKBB_coloc/PANUKBB_AMR_UKBB_AMR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/CSA_parquets --dir2 EST_UK_META/UKBB_CSA_parquets --results PANUKBB_coloc/PANUKBB_CSA_UKBB_CSA_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/EAS_parquets --dir2 EST_UK_META/UKBB_EAS_parquets --results PANUKBB_coloc/PANUKBB_EAS_UKBB_EAS_results.tsv --p12 1e-6 --H4 0.8"

  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/META_HQ_parquets\" --dir2 EST_UK_META/UKBB_EUR_parquets --results PANUKBB_coloc/PANUKBB_META_HQ_UKBB_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/META_HQ_parquets\" --dir2 EST_UK_META/EstBB_parquets --results PANUKBB_coloc/PANUKBB_META_HQ_EstBB_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/META_HQ_parquets\" --dir2 EST_UK_META/UKBB_AFR_parquets --results PANUKBB_coloc/PANUKBB_META_HQ_UKBB_AFR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/META_HQ_parquets\" --dir2 EST_UK_META/UKBB_MID_parquets --results PANUKBB_coloc/PANUKBB_META_HQ_UKBB_MID_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/META_HQ_parquets\" --dir2 EST_UK_META/UKBB_AMR_parquets --results PANUKBB_coloc/PANUKBB_META_HQ_UKBB_AMR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/META_HQ_parquets\" --dir2 EST_UK_META/UKBB_CSA_parquets --results PANUKBB_coloc/PANUKBB_META_HQ_UKBB_CSA_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/META_HQ_parquets\" --dir2 EST_UK_META/UKBB_EAS_parquets --results PANUKBB_coloc/PANUKBB_META_HQ_UKBB_EAS_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/META_HQ_parquets\" --dir2 EST_UK_META/meta_EUR_parquets --results PANUKBB_coloc/PANUKBB_META_HQ_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"

  "gpu-coloc -r --dir1 \"/path/to/PANUKBB_parquets_fixed/EUR_parquets\" --dir2 EST_UK_META/meta_EUR_parquets --results PANUKBB_coloc/PANUKBB_EUR_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/AFR_parquets --dir2 EST_UK_META/meta_EUR_parquets --results PANUKBB_coloc/PANUKBB_AFR_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/MID_parquets --dir2 EST_UK_META/meta_EUR_parquets --results PANUKBB_coloc/PANUKBB_MID_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/AMR_parquets --dir2 EST_UK_META/meta_EUR_parquets --results PANUKBB_coloc/PANUKBB_AMR_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/CSA_parquets --dir2 EST_UK_META/meta_EUR_parquets --results PANUKBB_coloc/PANUKBB_CSA_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
  "gpu-coloc -r --dir1 /path/to/PANUKBB_parquets_fixed/EAS_parquets --dir2 EST_UK_META/meta_EUR_parquets --results PANUKBB_coloc/PANUKBB_EAS_meta_EUR_results.tsv --p12 1e-6 --H4 0.8"
)

eval "${commands[$SLURM_ARRAY_TASK_ID]}"
