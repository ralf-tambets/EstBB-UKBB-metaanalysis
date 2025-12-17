#!/bin/bash

#SBATCH --time=180:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30GB
#SBATCH --job-name="MVP significant variants"
#SBATCH --partition=amd
#SBATCH --mail-user=mihkelje@ut.ee
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --array=0-3

module load python/3.12.3
source path/to/venv/bin/activate

datasets=(EUR META AMR AFR)

current_dataset=${datasets[$SLURM_ARRAY_TASK_ID]}

python3 MVP_significant_variants.py --pop ${current_dataset}

