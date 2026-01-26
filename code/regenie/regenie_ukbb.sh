#step 1
regenie \
--step 1 \
--bed $bfile_name \
--phenoFile $tsv_name \
--covarFile $tsv_name \
--extract $snplist_file_name \
--phenoColList $pheno_list \
--covarColList $covar_list \
--bsize 1000 \
--out regenie_step1 \
--qt \
--apply-rint \
--threads $cpus \
--gz \
--lowmem

#step 2
regenie \
--step 2 \
--bgen $bgen_file \
--sample $sample_file \
--ref-first \
--phenoFile $tsv_file \
--covarFile $tsv_file \
--phenoColList $pheno_list \
--covarColList $covar_list \
--qt \
--pred $pred_list \
--bsize 1000 \
--minMAC 20 \
--minINFO 0.6 \
--apply-rint \
--threads $cpus \
--out assoc_chr1 \
--gz