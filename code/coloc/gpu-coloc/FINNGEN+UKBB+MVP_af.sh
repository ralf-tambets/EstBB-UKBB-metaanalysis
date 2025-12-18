#!/bin/bash

#SBATCH --time=180:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30GB
#SBATCH --job-name="FMU get af"
#SBATCH --partition=amd
#SBATCH --mail-user=mihkelje@ut.ee
#SBATCH --mail-type=BEGIN,END,FAIL

module load python/3.12.3
source /path/to/venv/bin/activate # requires pip install gpu-coloc

python3 FINNGEN+UKBB+MVP_af.py

