#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(stringr)
library(susieR)
library(dotgen)
library(readr)
library(Rfast)
library(dplyr)
library(arrow)

#written by Ida Rahu
compute_yty <- function(beta, se, p, R, n, k) {
  # beta and se should be standarised
  beta_s <- beta * sqrt(2 * p * (1 - p))
  se_s <- se * sqrt(2 * p * (1 - p))
  
  # Y'Y =  Bj^2 (Xj'Xj) + Var(Bj)(Xj'Xj)(N - k)
  XjtXj <- (n - 1) * diag(R)
  yty <- beta_s**2 * XjtXj + se_s**2 * XjtXj * (n - k)
  
  return(median(yty))
}

#written by Ida Rahu
summarize.susie.cs <- function(object, orig_vars, R, ..., low_purity_threshold = 0.5) {
  if (is.null(object$sets)) {
    stop("Cannot summarize SuSiE object because credible set information is not available")
  }
  variables <- data.frame(cbind(1:length(object$pip), object$pip, -1, NA, NA, NA))
  colnames(variables) <- c("variable", "variable_prob", "cs", "cs_specific_prob", "low_purity", "lead_r2")
  rownames(variables) <- NULL
  added_vars <- c()
  if (object$null_index > 0) variables <- variables[-object$null_index, ]
  if (!is.null(object$sets$cs)) {
    cs <- data.frame(matrix(NA, length(object$sets$cs), 5))
    colnames(cs) <- c("cs", "cs_log10bf", "cs_avg_r2", "cs_min_r2", "variable")
    for (i in 1:length(object$sets$cs)) {
      if (any(object$sets$cs[[i]] %in% added_vars)) {
        print(
          sprintf("Skipping cs %d as there is an overlap between variants in this cs and previous credible sets", i)
        )
        print("Removed cs variants:")
        print(orig_vars[object$sets$cs[[i]], ], max = length(object$sets$cs[[i]]))
        next
      } else {
        added_vars <- append(added_vars, object$sets$cs[[i]])
      }
      in_cs_idx <- which(variables$variable %in% object$sets$cs[[i]])
      variables$cs[in_cs_idx] <- object$sets$cs_index[[i]]
      variables[in_cs_idx, "cs_specific_prob"] <- object$alpha[object$sets$cs_index[[i]], object$sets$cs[[i]]]
      variables$low_purity[in_cs_idx] <- object$sets$purity$min.abs.corr[i] < low_purity_threshold
      lead_pip_idx <- in_cs_idx[which.max(variables$variable_prob[in_cs_idx])]
      variables$lead_r2 <- R[lead_pip_idx, ]^2
      
      cs$cs[i] <- object$sets$cs_index[[i]]
      cs$cs_log10bf[i] <- log10(exp(object$lbf[cs$cs[i]]))
      cs$cs_avg_r2[i] <- object$sets$purity$mean.abs.corr[i]^2
      cs$cs_min_r2[i] <- object$sets$purity$min.abs.corr[i]^2
      cs$low_purity[i] <- object$sets$purity$min.abs.corr[i] < low_purity_threshold
      cs$variable[i] <- paste(object$sets$cs[[i]], collapse = ",")
    }
    variables <- variables[order(variables$variable_prob, decreasing = T), ]
  } else {
    cs <- NULL
  }
  return(list(vars = variables, cs = na.omit(cs)))
}

#based on work by Ida Rahu
susie_ss_wrapper <-function(df, R, n, L, estimate_residual_variance = TRUE, var_y = 1, prior_weights = NULL, min_abs_corr = 0.0, low_purity_threshold = 0.5) {
  beta <- df$beta
  se <- df$se
  fitted_bhat <- susie_rss(
    bhat = beta,
    shat = se,
    R = R,
    n = n,
    var_y = var_y,
    L = L,
    prior_weights = prior_weights,
    scaled_prior_variance = 0.1,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_variance = TRUE,
    standardize = TRUE,
    check_input = FALSE,
    min_abs_corr = min_abs_corr
  )
  
  cs_summary <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
  variables <-  as.data.frame(cs_summary$vars) %>%
    dplyr::rename(prob = variable_prob) %>%
    arrange(variable) %>%
    mutate(
      mean = susie_get_posterior_mean(fitted_bhat),
      sd = susie_get_posterior_sd(fitted_bhat)
    )
  cs <- cs_summary$cs
  
  cs_summary <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
  sets_95 <- fitted_bhat$sets
  fitted_bhat$sets <- susieR::susie_get_cs(fitted_bhat, coverage = 0.99, Xcorr = R, min_abs_corr = min_abs_corr)
  cs_summary_99 <- summarize.susie.cs(fitted_bhat, df, R, low_purity_threshold = low_purity_threshold)
  fitted_bhat$sets_99 <- fitted_bhat$sets
  fitted_bhat$sets <- sets_95
  
  variables_99 <-
    cs_summary_99$vars %>%
    dplyr::rename(prob = variable_prob) %>%
    arrange(variable) %>%
    mutate(
      mean = susie_get_posterior_mean(fitted_bhat),
      sd = susie_get_posterior_sd(fitted_bhat)
    )
  
  colnames(variables_99) <- paste0(colnames(variables_99),"_99")
  
  return(list(
    susie_obj = fitted_bhat,
    variables = variables,
    variables_99 = variables_99,
    cs = cs,
    cs_99 = cs_summary_99$cs
  ))
}

locus_chr = as.numeric(args[1])
locus_start = as.numeric(args[2])
locus_end = as.numeric(args[3])
LD_folder = args[4]
coloc_results_folder = args[5]
phenos=args[6]
n_samples=as.numeric(args[7])
gwas_folder=args[8]

#=======NB! MAY NEED TO CHANGE THESE=======
max_causal_SNPs = 10
n_covariates = 13
min_cs_corr = 0.5
low_purity_threshold = 0.5
GRCh = 38
gwas_suffix = "_EstBB.parquet"
LD_suffix = ".phased.vcor1"
#=======================================



print(paste(locus_chr, locus_start, locus_end))
print(paste(LD_folder, coloc_results_folder, gwas_folder))
print(phenos)
print(n_samples)
phenos <- unlist(strsplit(phenos, ";"))
message(length(phenos))

LD_file_prefix <- paste0(LD_folder, "/", locus_chr, "_", format(locus_start, scientific=FALSE), "_", format(locus_end, scientific=FALSE), LD_suffix)
LD_file <- paste0(LD_file_prefix, ".gz")
vars_file <- paste0(LD_file_prefix, ".vars")

message(LD_file)
message(vars_file)

if(!file.exists(LD_file)){
  message("no file")
  quit(save="no")
}

vars <- read_tsv(vars_file, col_names = F) %>% pull(1)
LD_matrix <- read.table(LD_file, fill = T, stringsAsFactors = F, col.names = vars, row.names = vars, check.names = F) #triangular LD matrix
dim(LD_matrix)

for (pheno in phenos){
  message(pheno)
  message(date())
  
  #expects parquet file NB! your suffix is probably different
  gwas_file <- paste0(gwas_folder, "/", pheno, "/", pheno, gwas_suffix)
  
  output_prefix <- paste0(coloc_results_folder, "/", pheno, "_", locus_chr, "_", format(locus_start, scientific=FALSE), "_", format(locus_end, scientific=FALSE))
  region <- paste0("chr", locus_chr,  ":", format(locus_start, scientific=FALSE), "-", format(locus_end, scientific=FALSE))
  
  if(file.exists(paste0(output_prefix, "_coloc5.tsv")) || file.exists(paste0(output_prefix, "_NULL_coloc5.tsv"))) {
    next
  }
  
  #expects parquet file
  gwas <- arrow::open_dataset(gwas_file) %>%
    filter(CHROM == locus_chr & between(GENPOS, locus_start, locus_end)) %>%
    filter(between(A1FREQ, 0.001, 0.999)) %>% #NB! MAF filter!
    collect() %>%
    arrange(CHROM, GENPOS)
  
  print(nrow(gwas))
  
  gwas <- gwas %>% #keep only variants that are in the LD matrix 
    filter(ID %in% vars) %>%
    mutate(maf = pmin(A1FREQ, 1-A1FREQ)) %>%
    rename(beta = BETA, se = SE)
  print(nrow(gwas))
  
  R <- as.matrix(LD_matrix[gwas$ID, gwas$ID]) #keep only variants that are in gwas 
  R[upper.tri(R)] <- t(R)[upper.tri(t(R))] #triangular to square matrix
  R[is.na(R)] <- 0 #replace na
  
  n <- n_samples
  L <- max_causal_SNPs
  prior_weights <- NULL
  
  gwas <- gwas %>%
    rename(rsid = ID, position = GENPOS, chromosome = CHROM, allele1 = ALLELE1, allele2 = ALLELE0)
  
  message('Step 1.')
  message(paste(length(colnames(R)), 'in analysis.'))
  
  
  message('Step 2.')
  yty <- compute_yty(beta = gwas$beta, se = gwas$se, p = gwas$maf, R = R, n = n, k = n_covariates)
  var_y <- yty / (n - 1)
  
  
  message('Step 3.')
  
  #added tryCatch - if it errors out, create empty files with _NULL_ in title
  res <- tryCatch(
    { 
      susie_ss_wrapper(df = gwas, R = R, n = n, L = L, estimate_residual_variance = TRUE, var_y = var_y, prior_weights = prior_weights, min_abs_corr = min_cs_corr, low_purity_threshold = low_purity_threshold)
    },
    error = function(cond) {
      message("ERROR")
      message("Here's the original error message:")
      message(conditionMessage(cond))
      # Choose a return value in case of error
      NULL
    }
  )
  
  if(is.null(res)){
    file.create(paste0(output_prefix, '_NULL_coloc5.tsv'))
    file.create(paste0(output_prefix, '_NULL_coloc3.tsv'))
    file.create(paste0(output_prefix, '_NULL_clpp.tsv'))
    next
  }
  
  
  message('Step 4.')
  susie_obj <- res$susie_obj
  
  lbf_variables <- susie_obj$lbf_variable
  transposed_lbf <- t(lbf_variables)
  transposed_lbf <- data.frame(rsid=row.names(transposed_lbf), transposed_lbf)
  transposed_lbf$old_position <- gwas$position
  transposed_lbf$A1 <- gwas$allele1
  transposed_lbf$A2 <- gwas$allele2
  
  old_rsid <- gwas$rsid
  old_position <- gwas$position
  gwas <- gwas %>% dplyr::rename('SNP'='rsid', 'CHR'='chromosome', 'BP'='position', 'A1'='allele1', 'A2'='allele2')
  gwas$rsid <- old_rsid
  gwas$old_position <- old_position
  gwas$CHR <- ifelse(gwas$CHR == 23, as.character('X'), gwas$CHR)
  
  if (GRCh == 38) {
    df_hg38 <- gwas
  } else {
    df_hg38 <- MungeSumstats::liftover(gwas, convert_ref_genome='hg38', ref_genome='hg19')
  }
  
  message(paste('Number of variants before liftover:', nrow(gwas)))
  message(paste('Number of variants after liftover:', nrow(df_hg38)))
  
  df_hg38$molecular_trait_id <- pheno
  df_hg38$region <- paste0('chr', unique(df_hg38$CHR),':', min(df_hg38$BP),'-', max(df_hg38$BP))
  message(paste('Old region:', region))
  message(paste('New region:', unique(df_hg38$region)))
  
  df_hg38$SNP <- paste0('chr', df_hg38$CHR, '_', df_hg38$BP, '_', df_hg38$A1, '_', df_hg38$A2)
  
  final_lbf <- inner_join(df_hg38, transposed_lbf, by=c('rsid', 'old_position', 'A1', 'A2'))
  final_lbf <- final_lbf[order(final_lbf$BP),]
  
  
  data_coloc5 <- final_lbf[, c('molecular_trait_id', 'region', 'SNP', 'CHR', 'BP', 'X1', 'X2',
                               'X3', 'X4', 'X5', 'X6', 'X7', 'X8', 'X9', 'X10')]
  
  colnames(data_coloc5) <- c('molecular_trait_id', 'region', 'variant', 'chromosome', 'position',
                             'lbf_variable1', 'lbf_variable2', 'lbf_variable3', 'lbf_variable4', 
                             'lbf_variable5', 'lbf_variable6', 'lbf_variable7', 'lbf_variable8', 
                             'lbf_variable9', 'lbf_variable10')
  
  write.table(data_coloc5, paste0(output_prefix, '_coloc5.tsv'), 
              row.names=F, sep='\t', quote=F)  
  message(paste('Coloc5 data is written.')) 
  
  data_coloc3 <- df_hg38[, c('molecular_trait_id', 'region', 'SNP', 'A1', 'A2', 'CHR', 'BP', 'maf',
                             'beta', 'se', 'LOG10P', 'INFO')]
  
  colnames(data_coloc3) <- c('molecular_trait_id', 'region', 'variant', 'ref', 'alt', 'chromosome',
                             'position', 'maf', 'beta', 'se', 'log10p', 'info')
  
  write.table(data_coloc3, paste0(output_prefix, '_coloc3.tsv'), 
              row.names=F, sep='\t', quote=F)  
  message(paste('Coloc3 data is written.'))
  
  df_rsid <- gwas$rsid
  variables <- cbind(df_rsid, res$variables)
  colnames(variables)[1] <- 'rsid'
  variables$old_position <- gwas$old_position
  variables$A1 <- gwas$A1
  variables$A2 <- gwas$A2
  variables$LOG10P <- gwas$LOG10P
  
  
  pip <- as.data.frame(t(susie_obj$alpha))
  colnames(pip) <- paste0("alpha", 1:10) #added alpha matrix to output - this are the exact PIP values
  pip$rsid <- rownames(pip)
  pip$pip <- susie_obj$pip
  pip$old_position <- gwas$old_position
  pip$A1 <- gwas$A1
  pip$A2 <- gwas$A2
  
  variables <- merge(variables, pip, by=c('rsid', 'old_position', 'A1', 'A2'))
  
  variables_clpp <- variables %>% filter(cs > 0)
  
  cs <- res$cs
  if(is.null(cs)) {
    write.table(variables_clpp, paste0(output_prefix, '_NULL_clpp.tsv'),
                row.names=F, sep='\t')
  } else {
    variables_clpp <- merge(variables_clpp, cs, by='cs')
    variables_clpp_hg38 <- inner_join(df_hg38, variables_clpp, by=c('rsid', 'old_position', 'A1', 'A2', "LOG10P"))
    
    message(paste('Number of variants in credible sets before liftover:', nrow(variables_clpp)))
    message(paste('Number of variants in credible sets after liftover:', nrow(variables_clpp_hg38)))
    
    variables_clpp_hg38 <- variables_clpp_hg38[order(variables_clpp_hg38$BP),]
    
    variables_clpp_hg38$z <- zsc(10^-(variables_clpp_hg38$LOG10P), variables_clpp_hg38$beta)
    variables_clpp_hg38$cs_index <- paste0('L', variables_clpp_hg38$cs)
    variables_clpp_hg38$cs_id <- paste0(variables_clpp_hg38$molecular_trait_id,
                                        '_', variables_clpp_hg38$region,
                                        '_', variables_clpp_hg38$cs_index)
    
    data_clpp <- variables_clpp_hg38[, c('molecular_trait_id', 'region', 'SNP', 'CHR', 'BP', 'A1', 'A2',
                                         'cs_id', 'cs_index', paste0("alpha", 1:10), 'pip', 'z')]
    
    colnames(data_clpp) <- c('molecular_trait_id', 'region', 'variant', 'chromosome', 'position',
                             'ref', 'alt', 'cs_id', 'cs_index', paste0("alpha", 1:10), 'pip', 'z')
    
    write.table(data_clpp, paste0(output_prefix, '_clpp.tsv'), 
                row.names=F, sep='\t', quote=F)
    
    message(paste('CLPP data is written.'))
    rm(R)
    gc()
  }
}

message(date())
message("Done")
