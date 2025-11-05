#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(stringr)
library(susieR)
library(dotgen)
library(readr)
library(dplyr)
library(arrow)
library(rbcor)

options(scipen=999)
print(date())

for(i in 1:13){
  message(args[i])
}

locus_chr = as.numeric(args[1])
lead_position = as.numeric(args[2])
locus_start = as.numeric(args[3])
locus_end = as.numeric(args[4])
pheno = args[5]
gwas_folder = args[6]
bcor_folder = args[7]
z_folder = args[8]
coloc_results_folder = args[9]
bcor_prefix = args[10]
gwas_suffix = args[11]
n_covariates = as.numeric(args[12])
is_metaanalysis = as.logical(args[13])
LD_window_size = as.numeric(args[14])

#=======NB! MAY NEED TO CHANGE THESE=======
prior_weights = NULL
independent_lead_LD_threshold = 0.05
max_causal_SNPs = 10
min_cs_corr = 0.5
low_purity_threshold = 0.5
maf_limit = 0.001
#=======================================

print(paste("Locus chr:", locus_chr))
print(paste("Lead position:", lead_position))
print(paste("Locus start:", locus_start))
print(paste("Locus end:", locus_end))
print(paste("Phenotype:", pheno))
print(paste("GWAS folder:", gwas_folder))
print(paste("BCOR folder:", bcor_folder))
print(paste("Z folder:", z_folder))
print(paste("Coloc results folder:", coloc_results_folder))
print(paste("BCOR prefix:", bcor_prefix))
print(paste("GWAS suffix:", gwas_suffix))
print(paste("Number of covariates:", n_covariates))
print(paste("Is metaanalysis:", is_metaanalysis))
print(paste("Prior weights:", prior_weights))
print(paste("Independent lead LD threshold:", independent_lead_LD_threshold))
print(paste("Max causal SNPs:", max_causal_SNPs))
print(paste("Min cs corr:", min_cs_corr))
print(paste("Low purity threshold:", low_purity_threshold))
print(paste("LD window size:", LD_window_size))
print(paste("MAF limit:", maf_limit))


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
susie_ss_wrapper <-function(df, R, n, L, estimate_residual_variance = TRUE, var_y = 1, prior_weights = NULL, min_abs_corr = 0.0, low_purity_threshold = 0.5, check_input = FALSE) {
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
    check_input = check_input,
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


finemap_start <- max(1, lead_position - LD_window_size)
finemap_end <- lead_position + LD_window_size

bcor_file_name <- paste0(bcor_prefix, "_chr", locus_chr, "_", format(locus_start, scientific = FALSE), "_", format(locus_end, scientific = FALSE))

z <- read_delim(paste0(z_folder, bcor_file_name, ".z")) %>%
  mutate(chromosome = as.numeric(gsub("chr", "", chromosome))) %>%
  mutate(variant = paste(chromosome, position, allele1, allele2, sep = "_")) 

if (is_metaanalysis){
  pheno_data <- arrow::open_dataset(paste0(gwas_folder, pheno, gwas_suffix)) %>%
    filter(CHROM == as.numeric(locus_chr)) %>%
    filter(between(GENPOS, finemap_start, finemap_end)) %>%
    mutate(A1FREQ = AC/(2*N)) %>%
    mutate(MAF = pmin(A1FREQ, 1-A1FREQ)) %>%
    filter(MAF > maf_limit) %>%
    select(chromosome = CHROM, position = GENPOS, allele1 = ALLELE0, allele2 = ALLELE1, N, beta = Effect, se = StdErr, LOG10P, maf = MAF, starts_with("INFO")) %>%
    collect() 
} else {
  pheno_data <- arrow::open_dataset(paste0(gwas_folder, pheno, "/", pheno, gwas_suffix)) %>%
    filter(CHROM == as.numeric(locus_chr)) %>%
    filter(between(GENPOS, finemap_start, finemap_end)) %>%
    mutate(MAF = pmin(A1FREQ, 1-A1FREQ)) %>%
    filter(MAF > maf_limit) %>%
    select(chromosome = CHROM, position = GENPOS, allele1 = ALLELE0, allele2 = ALLELE1, N, beta = BETA, se = SE, LOG10P, maf = MAF, INFO) %>%
    collect() 
}

z <- z %>%
  left_join(pheno_data, 
            by = join_by(chromosome, position, allele1, allele2),
            relationship = "one-to-one")

rm(pheno_data)

z_variants_to_keep <- which(between(z$position, lead_position - LD_window_size, lead_position + LD_window_size) & !is.na(z$beta))
z <- z[z_variants_to_keep,]
z <- z %>%
  mutate(molecular_trait_id = pheno,
         region = paste0('chr', ":", locus_chr, min(position), "-", max(position))) 

ld <- read_bcor(paste0("bcor_files/", bcor_file_name, ".bcor"))$read_corr()
ld <- ld[z_variants_to_keep, z_variants_to_keep]
colnames(ld) <- z$variant
rownames(ld) <- z$variant

gc()

output_prefix <- paste0(coloc_results_folder, "/", bcor_prefix, "_", pheno, "_", locus_chr, "_", finemap_start, "_", finemap_end)

n <- max(z$N)#NB
yty <- compute_yty(beta = z$beta, se = z$se, p = z$maf, R = ld, n = z$N, k = n_covariates)
var_y <- yty / (n - 1)


L <- max_causal_SNPs

#added tryCatch - if it errors out, create empty files with _NULL_ in title
res <- tryCatch(
  { 
    susie_ss_wrapper(df = z, 
                     R = ld, 
                     n = n, 
                     L = L, 
                     estimate_residual_variance = FALSE, #https://github.com/stephenslab/susieR/issues/162
                     var_y = var_y, 
                     prior_weights = prior_weights, 
                     min_abs_corr = min_cs_corr, 
                     low_purity_threshold = low_purity_threshold, 
                     check_input = FALSE) #kui true, siis XtX is not a positive semidefinite matrix. kui false, siis missing value where TRUE/FALSE needed (diagonaali sqrt läheb negatiivseks) Error in if (neg.loglik.logscale(lV, betahat = betahat, shat2 = shat2,
  },
  error = function(cond) {
    message(output_prefix)
    message(conditionMessage(cond))
    message("===")
    # Choose a return value in case of error
    NULL
  }
)

if(is.null(res)){
  file.create(paste0(output_prefix, '_NULL_coloc5.tsv'))
  file.create(paste0(output_prefix, '_NULL_coloc3.tsv'))
  file.create(paste0(output_prefix, '_NULL_clpp.tsv'))
  quit()
}

saveRDS(object = res, file = paste0(output_prefix, ".rds"))

data_coloc5 <- z %>%
  select(molecular_trait_id, region, variant, chromosome, position, ref = allele1, alt = allele2) %>%
  cbind(res$susie_obj$lbf_variable %>% t() %>% as.data.frame()) %>%
  arrange(chromosome, position, ref, alt) %>%
  rename_with(~ gsub("V", "lbf_variable", .x), starts_with("V"))

write_tsv(data_coloc5, paste0(output_prefix, '_coloc5.tsv'))

message(paste('Coloc5 data is written.')) 

data_coloc3 <- z %>%
  select(molecular_trait_id, region, variant, ref = allele1, alt = allele2, chromosome, position, maf, beta, se, LOG10P, starts_with("INFO"))

write_tsv(data_coloc3, paste0(output_prefix, '_coloc3.tsv'))

message(paste('Coloc3 data is written.'))

data_clpp <- z %>% 
  select(molecular_trait_id, region, variant, chromosome, position, beta, ref = allele1, alt = allele2, LOG10P) %>%
  cbind(res$susie_obj$alpha %>% t() %>% as.data.frame()) %>%
  rename_with(~ gsub("V", "alpha", .x), starts_with("V")) %>%
  mutate(pip = res$susie_obj$pip) %>%
  cbind(res$variables) %>%
  filter(cs > 0)

cs <- res$cs
if(is.null(cs)) {
  write_tsv(data_clpp, paste0(output_prefix, '_NULL_clpp.tsv'))
} else {
  data_clpp <- merge(data_clpp, cs, by='cs') %>%
    mutate(z = zsc(10^(-LOG10P), beta)) %>%
    mutate(cs_index = paste0("L", cs)) %>%
    mutate(cs_id = paste(molecular_trait_id, region, cs_index, sep = "_")) %>%
    select(molecular_trait_id, region, variant, chromosome, position, ref, alt, cs_id, cs_index, starts_with("alpha"), pip, z)
  
  write_tsv(data_clpp, paste0(output_prefix, '_clpp.tsv'))
  message(paste('CLPP data is written.'))
}

print(date())
