#!/bin/bash

#SBATCH --time=150:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30GB
#SBATCH --job-name="meta variants"
#SBATCH --partition=amd
#SBATCH --mail-user=mihkelje@ut.ee
#SBATCH --mail-type=BEGIN,END,FAIL


module load python/3.12.3
source path/to/venv/bin/activate

python3 get_meta_EUR_lead_metabolite_betas.py
