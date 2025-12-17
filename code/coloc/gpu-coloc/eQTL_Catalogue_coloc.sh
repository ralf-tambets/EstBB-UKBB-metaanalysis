#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH --job-name="Meta eQTL Catalogue coloc"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --array=0-7

datasets=(EstBB UKBB_AFR UKBB_AMR UKBB_CSA UKBB_EAS UKBB_EUR UKBB_MID meta_EUR)

current_dataset=${datasets[$SLURM_ARRAY_TASK_ID]}

module load python/3.12.3
source /path/to/venv/bin/activate
mkdir eQTL_catalogue_results

gpu-coloc -r --dir1 /path/to/eQTL_Catalogue/txrev_parquets --dir2 EST_UK_META/${current_dataset}_parquets --results eQTL_catalogue_results/txrev_${current_dataset}_results.tsv --p12 1e-6 --H4 0.8
gpu-coloc -r --dir1 /path/to/eQTL_Catalogue/ge_microarray_aptamer_parquets --dir2 EST_UK_META/${current_dataset}_parquets --results eQTL_catalogue_results/ge_microarray_aptamer_${current_dataset}_results.tsv --p12 1e-6 --H4 0.8
gpu-coloc -r --dir1 /path/to/eQTL_Catalogue/exon_parquets --dir2 EST_UK_META/${current_dataset}_parquets --results eQTL_catalogue_results/exon_${current_dataset}_results.tsv --p12 1e-6 --H4 0.8
gpu-coloc -r --dir1 /path/to/eQTL_Catalogue/tx_parquets --dir2 EST_UK_META/${current_dataset}_parquets --results eQTL_catalogue_results/tx_${current_dataset}_results.tsv --p12 1e-6 --H4 0.8
gpu-coloc -r --dir1 /path/to/eQTL_Catalogue/leafcutter_parquets --dir2 EST_UK_META/${current_dataset}_parquets --results eQTL_catalogue_results/leafcutter_${current_dataset}_results.tsv --p12 1e-6 --H4 0.8
