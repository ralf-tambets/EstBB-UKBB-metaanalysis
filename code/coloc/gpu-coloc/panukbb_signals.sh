#!/bin/bash

#SBATCH --time=140:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="PANUKBB format"
#SBATCH --partition=amd
#SBATCH --array=0-6

module load python/3.12.3
source /path/to/venv/bin/activate

mkdir PANUKBB

datasets=(AFR AMR CSA EAS EUR MID META_HQ)

current_dataset=${datasets[$SLURM_ARRAY_TASK_ID]}

python3 panukbb_signals.py --pop ${current_dataset} --root /path/to/PAN_UKBB_lifted --out PANUKBB

gpu-coloc -f \
    --input          "PANUKBB/${current_dataset}_signals" \
    --output         "PANUKBB/${current_dataset}_parquets" \
    --input_summary  "PANUKBB/${current_dataset}_summary.tsv"