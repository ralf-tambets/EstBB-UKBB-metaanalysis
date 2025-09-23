#!/bin/bash

#SBATCH --time=140:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="format Est UKBB meta-analysis"
#SBATCH --partition=amd
#SBATCH --array=0-6

module load python/3.12.3
source /path/to/venv/bin/activate # requires pip install gpu-coloc
mkdir EST_UK_META

datasets=(EstBB UKBB_AFR UKBB_AMR UKBB_CSA UKBB_EAS UKBB_EUR UKBB_MID)

current_dataset=${datasets[$SLURM_ARRAY_TASK_ID]}

python3 meta_signals.py --pop ${current_dataset}

gpu-coloc -f \
    --input          "EST_UK_META/${current_dataset}_signals" \
    --output         "EST_UK_META/${current_dataset}_parquets" \
    --input_summary  "EST_UK_META/${current_dataset}_summary.tsv"