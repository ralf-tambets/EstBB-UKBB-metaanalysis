#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --mem=64GB
#SBATCH --job-name="LD"
#SBATCH --partition=main

module load plink2/devel

chr=$1
start=$2
end=$3
outname=${chr}_${start}_${end}_N1
pfile_prefix=...
keep_file=...

echo $chr $start $end
plink2 \
--pfile ${pfile_prefix}${chr} \
--keep $keep_file \
--rm-dup force-first \
--r-phased triangle ref-based \
--chr $chr \
--from-bp $start \
--to-bp $end \
--out $outname \
--memory 64000 \
--maf 0.001 \
--threads 16