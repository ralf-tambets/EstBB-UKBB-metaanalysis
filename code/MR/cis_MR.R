library(dplyr)
library(rbcor)
library(mrlocus)
library(MendelianRandomization)

#df - data frame of leads to LD prune, contains a column named "variant" that has values of underscore separated chromosome, position, ref, alt values (chromosome_position_ref_alt).
#bcor_file - bcor file to open
#variant_col - name of the column that contains a unique identifier for each variant
#if_r_matrix - whether the bcor file contains r values (if it does, this function squares them)
read_LD_from_bcor <- function(lead_df, bcor_file, if_r_matrix = TRUE) {
  bcor <- read_bcor(bcor_file)
  
  meta <- bcor$get_meta() %>%
    mutate(variant = paste(gsub("chr", "", chromosome), position, allele1, allele2, sep = "_")) 
  
  leads_idx <- which(meta$variant %in% lead_df$variant)
  
  ld <- bcor$read_corr(snps = leads_idx, snps2 = leads_idx)
  rownames(ld) <- meta$variant[leads_idx]
  colnames(ld) <- meta$variant[leads_idx]
  
  if(if_r_matrix){
    ld = ld^2
  }
  
  return(ld)
}

#leads_df - data frame of leads to LD prune
#ld_matrix - matrix of r2 values - should have row and column names corresponding to variant_col
#col_to_compare - name of the column in leads_df that is under review (e.g. pip or LOG10P)
#sort_descending - whether to pick the highest value in col_to_compare
#variant_col - name of the column that contains a unique identifier for each variant (should match with row and column names of ld_matrix)
#r2_limit - prune out all variants in higher LD than this
prune_by_LD <- function(lead_df, ld_matrix, col_to_compare = "pip", sort_descending = TRUE, variant_col = "variant", r2_limit = 0.01){
  sorted_df <- lead_df[order(lead_df[[col_to_compare]], decreasing = sort_descending),]
  
  pruned_df <- data.frame()
  
  if(!all(sorted_df[[variant_col]] %in% colnames(ld_matrix))) message("WARNING! Some lead variants are not present in the LD matrix!")
  
  while(nrow(sorted_df) > 0){
    lead_SNP = sorted_df[1,]
    
    snps_in_LD <- names(ld_matrix[lead_SNP[[variant_col]],][ld_matrix[lead_SNP[[variant_col]],] > r2_limit])
    
    pruned_df <- pruned_df %>%
      rbind(lead_SNP)
    
    sorted_df <- sorted_df[!(sorted_df[[variant_col]] %in% snps_in_LD),]
  }
  
  return(pruned_df)
}

#df - data frame with these columns: variant, beta_exposure, se_exposure, beta_outcome, se_outcome
#both betas should be multiplied with the sign of beta_exposure beforehand (to get all values to the right of the y-axis)
run_mendelianrandomization <- function(df, exposure_name = "exposure", outcome_name = "outcome"){
  if(nrow(df) < 2){
    message("Not enough lead variants")
    return(NULL)
  }
  
  MR_input <- MendelianRandomization::mr_input(
    exposure = exposure_name, 
    bx = df$beta_exposure, 
    bxse = df$se_exposure, 
    outcome = outcome_name,
    by = df$beta_outcome, 
    byse = df$se_outcome, 
    snps = df$variant
  )
  
  mr_ivw_results <- tryCatch({
    MendelianRandomization::mr_ivw(MR_input, alpha = 0.05, weights = "delta", model = "random")

  }, error = function(err) {
    print(paste("ERROR:  ",err))
    return(NULL)
  })
  
  return(mr_ivw_results)
}

#df - data frame with these columns: variant, beta_exposure, se_exposure, beta_outcome, se_outcome
#both betas should be multiplied with the sign of beta_exposure beforehand (to get all values to the right of the y-axis)
run_mrlocus <- function(df) {
  res <- list(beta_hat_a = df$beta_exposure,
              beta_hat_b = df$beta_outcome,
              sd_a=df$se_exposure,
              sd_b=df$se_outcome)
  
  res <- mrlocus::extractForSlope(res, plot = F)
  
  res <- mrlocus::fitSlope(res, iter=10000)
  
  return(res)
}

#res - output of run_mrlocus()
#main - title of the plot, usually the name of the gene/locus
#a - outcome name
#b - exposure name
#mrlocus::plotMrlocus(res, main, a, b)

#df - data frame with these columns: variant, beta_exposure, se_exposure, beta_outcome, se_outcome
#both betas should be multiplied with the sign of beta_exposure beforehand (to get all values to the right of the y-axis)
#res - output of run_mendelianrandomization()
#main - title of the plot, usually the name of the gene/locus
#xlab - outcome name
#ylab - exposure name
plot_mr_ivw_with_mrlocus <- function(df, res, locus_colors = NULL, q=c(.025,.975), sigma_mult=1.96, xlim=NULL, ylim=NULL, legend=TRUE, digits=3, col_slope="blue", col_band=rgb(0,0,1,.1), col_dashed=rgb(0,0,1,.5), legend_where = NULL, ...) {
  beta_hat_a = df$beta_exposure
  beta_hat_b = df$beta_outcome
  sd_a = df$se_exposure
  sd_b = df$se_outcome
  alpha.hat = res@Estimate
  sigma.hat = res@StdError
  
  stopifnot(length(q) == 2)
  qs <- paste0(q * 100,"%")
  signs = sign(beta_hat_a)
  beta_hat_a = beta_hat_a * signs
  beta_hat_b = beta_hat_b * signs
  if (is.null(xlim)) {
    xx <- max(beta_hat_a)
    xlim <- c(0, 1.5*xx)
  } else {
    xx <- 1.33*xlim[2]
  }
  yy <- 1.5*max(abs(beta_hat_b))
  if (is.null(ylim)) {
    ylim <- c(-yy, yy)
  }
  plot(beta_hat_a, beta_hat_b,
       xlim=xlim, ylim=ylim, type="n", ...)
  
  alpha.qs <- c(alpha.hat - sigma_mult*sigma.hat, alpha.hat + sigma_mult*sigma.hat)
  
  segments(0, 0, 2*xx, alpha.hat*2*xx, col=col_slope, lwd=2)
  segments(0, 0, 2*xx, (alpha.qs[1])*2*xx,
           col=col_dashed, lwd=2, lty=2)
  segments(0, 0, 2*xx, (alpha.qs[2])*2*xx,
           col=col_dashed, lwd=2, lty=2)
  abline(h=0, col=rgb(0,0,0,.25))
  
  arrows(beta_hat_a - sd_a, beta_hat_b,
         beta_hat_a + sd_a, beta_hat_b,
         code=3, angle=90, length=.05)
  arrows(beta_hat_a, beta_hat_b - sd_b,
         beta_hat_a, beta_hat_b + sd_b,
         code=3, angle=90, length=.05)
  
  if (is.null(locus_colors)){
    points(beta_hat_a, beta_hat_b, pch=19)
  } else {
    points(beta_hat_a, beta_hat_b, pch=19, col = locus_colors)
  }
  
  if (legend) {
    if (is.null(legend_where)){
      where <- if (alpha.hat > 0) "bottomleft" else "topleft"
    } else {
      where <- legend_where
    }
    
    slope.leg <- as.expression(bquote(paste("Estimate"," = ",
                                            .(round(alpha.hat,digits)))))
    slope.int.leg <- paste0(100*diff(q), "% CI = (",
                            round(alpha.qs[1],digits),
                            ", ",round(alpha.qs[2],digits),
                            ")")
    sigma.leg <- as.expression(bquote(paste("StdError"," = ",
                                            .(round(sigma.hat,digits)))))
    
    legend(where,
           lwd=c(2,2,5),
           lty=c(1,2,1),
           col=c(col_slope, col_dashed, col_band),
           inset=.05,
           y.intersp=1.1,
           bg="white",
           legend=c(slope.leg, slope.int.leg, sigma.leg))
  }
}

#example
#plot_mr_ivw_with_mrlocus(df = ..., res = ..., main = "Title", xlab = "Effect size of exposure", ylab = "Effect size of outcome")