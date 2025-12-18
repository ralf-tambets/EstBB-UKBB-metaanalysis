#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=60GB
#SBATCH --job-name="GPU Coloc Suzuki Aragam signals"
#SBATCH --partition=amd

module load python/3.12.3
source gpu-coloc/coloc_env/bin/activate

python3 Suzuki_Aragam_signals.py --summary Suzuki_Aragam_summary.tsv --output Suzuki_Aragam_signals --input_table diaseas_gwas_table.tsv

gpu-coloc -f \
    --input          "Suzuki_Aragam_signals" \
    --output         "Suzuki_Aragam_parquets" \
    --input_summary  "Suzuki_Aragam_summary.tsv"