#step 1
regenie \
    --step 1 \
    --bed $bed_prefix \
    --phenoFile $phenotype_file \
    --phenoColList $phenotype_id \
    --covarFile $covariate_file \
    --covarColList $covariate_list \
    --use-relative-path \
    --bsize 1000 \
    --lowmem \
    --apply-rint \
    --lowmem-prefix $prefix_tmp_rg \
    --gz \
    --threads $cpus \
    --out $prefix 

#step 2
regenie \
    --step 2 \
    --bgen $bgen \
    --sample $sample \
    --ref-first \
    --phenoFile $phenotype_file \
    --phenoColList $phenotype_id \
    --covarFile $covariate_file \
    --covarColList $covariate_list \
    --chr $chromosome \
    --minINFO 0.6 \
    --minMAC 20 \
    --bsize 1000 \
    --apply-rint \
    --threads $cpus \
    --pred $pred_list \
    --gz \
    --out ${prefix}_${chromosome}"\