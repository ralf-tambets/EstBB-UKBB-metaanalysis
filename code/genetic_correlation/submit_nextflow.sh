#!/bin/bash

#SBATCH --time=...
#SBATCH -N ...
#SBATCH --ntasks-per-node=...
#SBATCH --mem=...
#SBATCH --job-name="gc"

module load any/jdk/1.8.0_265
module load nextflow
module load any/singularity/3.5.3
module load squashfs/4.4

nextflow run main.nf -entry between_metabolites -resume