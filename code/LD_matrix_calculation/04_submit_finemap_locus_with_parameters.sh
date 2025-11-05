#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16GB
#SBATCH --job-name="finemap"
#SBATCH --partition=main

locus_chr=$1
lead_position=$2
locus_start=$3
locus_end=$4
pheno=$5
gwas_folder="/gpfs/helios/projects/alasoo/metabolite_gwas/EstBB_UKB_meta-analysis/UKBB_EUR/STEP2/"
bcor_folder="./bcor_files/"
z_folder="./z_files/"
coloc_results_folder="./susie_files/"
bcor_prefix="EUR"
gwas_suffix="_EUR.parquet"
n_covariates=23
is_metaanalysis="FALSE"
LD_window_size=1000000


module load any/jdk/1.8.0_265;
module load squashfs/4.4;
module load any/singularity/3.5.3;
echo $locus_chr $lead_position $locus_start $locus_end $pheno $gwas_folder $bcor_folder $z_folder $coloc_results_folder $bcor_prefix $gwas_suffix $n_covariates $is_metaanalysis $LD_window_size
singularity exec -B /gpfs/:/gpfs/ /gpfs/helios/home/rtambets/LDSTORE_DNANEXUS/quay.io-rtambets-susie-v25.11.2.img Rscript finemap_locus.R \
$locus_chr \
$lead_position \
$locus_start \
$locus_end \
$pheno \
$gwas_folder \
$bcor_folder \
$z_folder \
$coloc_results_folder \
$bcor_prefix \
$gwas_suffix \
$n_covariates \
$is_metaanalysis \
$LD_window_size;




