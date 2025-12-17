#!/bin/bash

#SBATCH --time=140:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12GB
#SBATCH --job-name="UKBPPP format"
#SBATCH --partition=amd
#SBATCH --array=0-6

module load python/3.12.3
source /path/to/venv/bin/activate

datasets=(
  "African"
  "Central_South_Asian"
  "East_Asian"
  "Middle_East"
  "American"
  "Combined"
  "European_(discovery)"
)

current_dataset=${datasets[$SLURM_ARRAY_TASK_ID]}

python3 UKBPPP_signals.py --pop ${current_dataset}

gpu-coloc -f \
    --input          "UKBPPP_gpu/${current_dataset}_signals" \
    --output         "UKBPPP_gpu/${current_dataset}_parquets" \
    --input_summary  "UKBPPP_gpu/${current_dataset}_summary.tsv"