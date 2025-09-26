#!/bin/bash

#SBATCH --time=140:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="MVP format"
#SBATCH --partition=amd
#SBATCH --array=0-3

module load python/3.12.3
source /path/to/venv/bin/activate # requires pip install gpu-coloc

mkdir MVP

datasets=(EUR META AMR AFR)

current_dataset=${datasets[$SLURM_ARRAY_TASK_ID]}

python3 MVP_signals.py --pop ${current_dataset}

gpu-coloc -f \
    --input          "MVP/${current_dataset}_signals" \
    --output         "MVP/${current_dataset}_parquets" \
    --input_summary  "MVP/${current_dataset}_summary.tsv"