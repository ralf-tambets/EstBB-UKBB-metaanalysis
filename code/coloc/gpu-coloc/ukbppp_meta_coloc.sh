#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="UKBPPP coloc, meta_EUR"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

#SBATCH --array=0-6

to_test=(
  "African"
  "Central_South_Asian"
  "East_Asian"
  "Middle_East"
  "American"
  "Combined"
  "European_(discovery)"
)

current_test_dataset=${to_test[$SLURM_ARRAY_TASK_ID]}

module load python/3.12.3
source /path/to/venv/bin/activate
mkdir UKBPPP_coloc
gpu-coloc -r --dir1 UKBPPP_gpu/${current_test_dataset}_parquets --dir2 EST_UK_META/meta_EUR_parquets --results UKBPPP_coloc/UKBPPP_${current_test_dataset}_meta_EUR_results.tsv --p12 1e-6 --H4 0
