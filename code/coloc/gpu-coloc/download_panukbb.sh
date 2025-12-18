#!/bin/bash

#SBATCH --time=150:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60GB
#SBATCH --job-name="Download PANUKBB"
#SBATCH --partition=amd



module load python/3.12.3
source /path/to/venv/bin/activate

mkdir /gpfs/space/projects/genomic_references/summary_stats/PAN_UKBB

python3 download_panukbb.py
