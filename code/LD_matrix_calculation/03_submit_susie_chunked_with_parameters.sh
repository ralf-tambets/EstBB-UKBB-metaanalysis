#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --mem=32GB
#SBATCH --job-name="fine-map"
#SBATCH --partition=main

module load any/jdk/1.8.0_265
module load squashfs/4.4
module load any/singularity/3.5.3

chr=$1
start=$2
end=$3
phenos=$4
LD_folder=...
results_folder=...
n_samples=...
gwas_folder=...

echo $chr $start $end $LD_folder $results_folder $phenos $n_samples $gwas_folder
singularity exec -B /gpfs/:/gpfs/ .../quay.io-rtambets-susie-v25.08.2.img Rscript 03_susie_chunked.R \
$chr \
$start \
$end \
$LD_folder \
$results_folder \
$phenos \
$n_samples \
$gwas_folder