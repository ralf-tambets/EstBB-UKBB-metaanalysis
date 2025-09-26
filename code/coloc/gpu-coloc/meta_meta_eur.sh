#!/bin/bash

#SBATCH --time=140:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="format Est UKBB meta-analysis meta_EUR"
#SBATCH --partition=amd


module load python/3.12.3
source /path/to/venv/bin/activate # requires pip install gpu-coloc

python3 meta_meta_eur.py

gpu-coloc -f \
    --input          "EST_UK_META/meta_EUR_signals" \
    --output         "EST_UK_META/meta_EUR_parquets" \
    --input_summary  "EST_UK_META/meta_EUR_wo_dupes.tsv"