#!/bin/bash

#SBATCH --time=180:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH --job-name="FinnGen+MVP+UKBB format"
#SBATCH --partition=amd

module load python/3.12.3
source /path/to/venv/bin/activate # requires pip install gpu-coloc

mkdir FinnGen+MVP+UKBB

python3 FINNGEN+UKBB+MVP.py

gpu-coloc -f \
    --input          "FinnGen+MVP+UKBB/FinnGen+MVP+UKBB_signals" \
    --output         "FinnGen+MVP+UKBB/FinnGen+MVP+UKBB_parquets" \
    --input_summary  "FinnGen+MVP+UKBB/FinnGen+MVP+UKBB_summary.tsv"